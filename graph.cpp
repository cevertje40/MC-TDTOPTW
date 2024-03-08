#include "graph.h"

void Node::addarc(Link* ARC)//dereference arc=>node
{
	ARC->next = first;
	first = ARC;
}

void Node::addprevious(Link* ARC)
{
	ARC->nextbackward = firstbackward;
	firstbackward = ARC;
}

void Nodep::addarc(Link* ARC)//dereference arc=>node
{
	ARC->next = first;
	first = ARC;
}

void Nodep::addprevious(Link* ARC)
{
	ARC->nextbackward = firstbackward;
	firstbackward = ARC;
}


Graph::Graph(int maxnodes, int maxlinks, int maxtimeslots) 
{
	n.resize(maxnodes);
	np.resize(maxnodes);
	ifstream file("..\\..\\datasets\\bemobile\\nodes_cleaned.csv");
	std::string value;
	int counter = 0;
	bool stop = false;
	while (stop == false)
	{
		getline(file, value, ';');
		n[counter].id = atoi(value.c_str());
		np[counter].id = atoi(value.c_str());
		getline(file, value, ';');
		char* pEnd;
		n[counter].longitude = strtod(value.c_str(), &pEnd);
		np[counter].longitude = strtod(value.c_str(), &pEnd);
		getline(file, value, '\n');
		n[counter].latitude = strtod(value.c_str(), &pEnd);
		np[counter].latitude = strtod(value.c_str(), &pEnd);
		++counter;
		if (file.peek() == EOF)
			stop = true;
	}
	file.close();
	/*************
	reading link data
	***************/
	l.resize(maxlinks);
	file.open("..\\..\\datasets\\bemobile\\links_cleaned.csv");
	counter = 0;
	stop = false;
	while (stop == false)
	{
		getline(file, value, ';');
		l[counter].link_id = atoi(value.c_str());
		getline(file, value, ';');
		l[counter].from = atoi(value.c_str()) - 1;
		getline(file, value, ';');
		l[counter].to = atoi(value.c_str()) - 1;
		getline(file, value, '\n');
		l[counter].optimaltt = double(atoi(value.c_str())) / (3600 * 1000);
		if (file.peek() == EOF)//laatste ; wegdoen
			stop = true;
		n[l[counter].from].addarc(&l[counter]);//ad forward link to startnode
		n[l[counter].to].addprevious(&l[counter]);//ad backward link to endnode
		np[l[counter].from].addarc(&l[counter]);//ad forward link to startnode
		np[l[counter].to].addprevious(&l[counter]);//ad backward link to endnode
		++counter;
	}//end second while
	file.close();
	/*************
	reading traveltime data
	***************/
	file.open("..\\..\\datasets\\bemobile\\traveltime_cleaned.csv");
	for (int i = 0; i < maxlinks; ++i)
	{
		l[i].traveltime.resize(maxtimeslots + 1);
		for (int t = 0; t < maxtimeslots + 1; ++t)
		{
			getline(file, value, '\n');
			l[i].traveltime[t] = double(atoi(value.c_str())) / (3600 * 1000);
		}
	}
	file.close();
	//create gamma link distributions
	#pragma omp parallel
	{//start parallel session
	#pragma omp for nowait
		for (int i = 0; i < maxlinks; ++i)
		{
			l[i].k.resize(maxtimeslots + 1);
			l[i].th.resize(maxtimeslots + 1);
			for (int t = 0; t < maxtimeslots + 1; ++t)
			{
				if (l[i].traveltime[t] - l[i].optimaltt > 0)
				{
					double sd = calculate_stdv_pop_d(calculate_mean_d(l[i].traveltime), l[i].traveltime);
					//mean is the estimate at time t here not the mean over all periods t
					l[i].k[t] = pow(l[i].traveltime[t], 2) / (pow(sd, 2));
					l[i].th[t] = pow(sd, 2) / l[i].traveltime[t];
				}
				else
				{//divide by a very small sd to have meaningfull numbers
					l[i].k[t] = pow(l[i].traveltime[t], 2) / (pow(0.01, 2));
					l[i].th[t] = pow(0.01, 2) / l[i].traveltime[t];
				}
			}//for all timeslots
		}// for all links
	}//end parallel
	cout << "be mobile roadnetwork data succesfully read in" << endl;
}


inline  int BinaryMinHeap::getLeftChildIndex(int nodeIndex)
{
	return nodeIndex << 1; //2 * nodeIndex;
}
inline   int BinaryMinHeap::getRightChildIndex(int nodeIndex)
{
	return 1 + (nodeIndex << 1); //2 * nodeIndex + 1;
}
inline   int BinaryMinHeap::getParentIndex(int nodeIndex)
{
	return nodeIndex >> 1; //(nodeIndex/ 2);
}


BinaryMinHeap::BinaryMinHeap(int size) //constructor
{
	Nodes = new Node * [size + 1];
	heapSize = 0;
}

BinaryMinHeap::~BinaryMinHeap() //destructor
{
	delete[] Nodes;
	Nodes = NULL;
}

void BinaryMinHeap::siftUp(int nodeIndex)
{
	Nodes[0] = Nodes[nodeIndex];
	int parentIndex = getParentIndex(nodeIndex);
	while (nodeIndex > 1 && Nodes[parentIndex]->score > Nodes[0]->score)
	{
		Nodes[nodeIndex] = Nodes[parentIndex];
		Nodes[nodeIndex]->hdex = nodeIndex;
		nodeIndex = parentIndex;
		parentIndex = getParentIndex(nodeIndex);

	}
	Nodes[nodeIndex] = Nodes[0];
	Nodes[nodeIndex]->hdex = nodeIndex;
}

void BinaryMinHeap::insert(Node* Node)
{
	heapSize++;
	Nodes[heapSize] = Node;
	Nodes[heapSize]->hdex = heapSize;
	siftUp(heapSize);
}

void BinaryMinHeap::siftDown(int nodeIndex)
{
	int leftChildIndex, rightChildIndex, minIndex;
	Node* temp;
	leftChildIndex = getLeftChildIndex(nodeIndex);
	rightChildIndex = getRightChildIndex(nodeIndex);
	if (rightChildIndex > heapSize)
	{
		if (leftChildIndex > heapSize)
			return;
		else
			minIndex = leftChildIndex;
	}
	else
	{
		if (Nodes[leftChildIndex]->score <= Nodes[rightChildIndex]->score)
			minIndex = leftChildIndex;
		else
			minIndex = rightChildIndex;
	}
	if (Nodes[nodeIndex]->score > Nodes[minIndex]->score)
	{
		temp = Nodes[minIndex];
		Nodes[minIndex] = Nodes[nodeIndex];
		Nodes[nodeIndex] = temp;
		Nodes[nodeIndex]->hdex = nodeIndex;
		Nodes[minIndex]->hdex = minIndex;
		siftDown(minIndex);
	}
}

Node* BinaryMinHeap::extractMin()
{
	Nodes[0] = Nodes[1];//store min to temp
	Nodes[1]->hdex = -1;
	Nodes[1] = Nodes[heapSize];
	Nodes[1]->hdex = 1;
	heapSize--;
	if (heapSize > 0)
		siftDown(1);
	return Nodes[0];
}//end extractmin


Node* BinaryMinHeap::extractMintwee()
{
	if (heapSize > 0)
	{
		Nodes[0] = Nodes[1];//store min to temp
		Nodes[1]->hdex = -1;
		Nodes[1] = Nodes[heapSize];
		Nodes[1]->hdex = 1;
		heapSize--;
		siftDown(1);
		return Nodes[0];
	}
	else
	{
		return NULL;
	}

}//end extractmintwee
