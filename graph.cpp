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

inline int BinaryMinHeap::getLeftChildIndex(int nodeIndex)
{
	return nodeIndex << 1; //2 * nodeIndex;
}
inline int BinaryMinHeap::getRightChildIndex(int nodeIndex)
{
	return 1 + (nodeIndex << 1); //2 * nodeIndex + 1;
}
inline int BinaryMinHeap::getParentIndex(int nodeIndex)
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

inline  int BinaryMinHeaps::getLeftChildIndex(int nodeIndex)
{
	return nodeIndex << 1; //2 * nodeIndex;
}

inline   int BinaryMinHeaps::getRightChildIndex(int nodeIndex)
{
	return 1 + (nodeIndex << 1); //2 * nodeIndex + 1;
}

inline   int BinaryMinHeaps::getParentIndex(int nodeIndex)
{
	return nodeIndex >> 1; //(nodeIndex/ 2);
}

BinaryMinHeaps::BinaryMinHeaps() // default constructor
{
}

BinaryMinHeaps::BinaryMinHeaps(int size, int thread) //constructor
{
	Nodes[thread] = new Nodep * [size + 1];
	heapSize[thread] = 0;
}

BinaryMinHeaps::~BinaryMinHeaps() //destructor
{
}

void BinaryMinHeaps::free(int thread)//threaded destructor
{
	delete[] Nodes[thread];
	Nodes[thread] = NULL;
}

void BinaryMinHeaps::siftUp(int nodeIndex, int thread)
{
	Nodes[thread][0] = Nodes[thread][nodeIndex];
	int parentIndex = getParentIndex(nodeIndex);

	while (nodeIndex > 1 && Nodes[thread][parentIndex]->score[thread] > Nodes[thread][0]->score[thread])
	{
		Nodes[thread][nodeIndex] = Nodes[thread][parentIndex];
		Nodes[thread][nodeIndex]->hdex[thread] = nodeIndex;
		nodeIndex = parentIndex;
		parentIndex = getParentIndex(nodeIndex);
	}
	Nodes[thread][nodeIndex] = Nodes[thread][0];
	Nodes[thread][nodeIndex]->hdex[thread] = nodeIndex;
}

void BinaryMinHeaps::insert(Nodep* Node, int thread)
{
	heapSize[thread]++;
	Nodes[thread][heapSize[thread]] = Node;
	Nodes[thread][heapSize[thread]]->hdex[thread] = heapSize[thread];
	siftUp(heapSize[thread], thread);
}

void BinaryMinHeaps::siftDown(int nodeIndex, int thread)
{
	int leftChildIndex, rightChildIndex, minIndex;
	Nodep* temp;
	leftChildIndex = getLeftChildIndex(nodeIndex);
	rightChildIndex = getRightChildIndex(nodeIndex);
	if (rightChildIndex > heapSize[thread])
	{
		if (leftChildIndex > heapSize[thread])
			return;
		else
			minIndex = leftChildIndex;
	}
	else
	{
		if (Nodes[thread][leftChildIndex]->score[thread] <= Nodes[thread][rightChildIndex]->score[thread])
			minIndex = leftChildIndex;
		else
			minIndex = rightChildIndex;
	}
	if (Nodes[thread][nodeIndex]->score[thread] > Nodes[thread][minIndex]->score[thread])
	{
		temp = Nodes[thread][minIndex];
		Nodes[thread][minIndex] = Nodes[thread][nodeIndex];
		Nodes[thread][nodeIndex] = temp;
		Nodes[thread][nodeIndex]->hdex[thread] = nodeIndex;
		Nodes[thread][minIndex]->hdex[thread] = minIndex;
		siftDown(minIndex, thread);
	}
}

Nodep* BinaryMinHeaps::extractMin(int thread)
{
	Nodes[thread][0] = Nodes[thread][1];//store min to temp
	Nodes[thread][1]->hdex[thread] = -1;
	//cout<<"final score for"<<node_id[index[1]]<<":"<<node_score[index[1]]<<endl;
	Nodes[thread][1] = Nodes[thread][heapSize[thread]];
	Nodes[thread][1]->hdex[thread] = 1;
	heapSize[thread]--;
	if (heapSize[thread] > 0)
		siftDown(1, thread);
	return Nodes[thread][0];
	//cout<<"check the heap"<<endl;
}//end extractmin

Nodep* BinaryMinHeaps::extractMintwee(int thread)
{
	if (heapSize[thread] > 0)
	{
		Nodes[thread][0] = Nodes[thread][1];//store min to temp
		Nodes[thread][1]->hdex[thread] = -1;
		Nodes[thread][1] = Nodes[thread][heapSize[thread]];
		Nodes[thread][1]->hdex[thread] = 1;
		heapSize[thread]--;
		siftDown(1, thread);
		return Nodes[thread][0];
	}
	else
	{
		return NULL;
	}

}//end extractmintwee

double Graph::calculate_mean_d(vector<double>& input)//mean
{
	double sum = 0;
	for (int i = 0; i < input.size(); i++)
	{
		sum += input[i];
	}
	return (sum / input.size());
}

double Graph::calculate_stdv_pop_d(double mean, vector<double>& input)//population standard deviation
{
	double temp = 0;
	for (int i = 0; i < input.size(); i++)
	{
		temp += (input[i] - mean) * (input[i] - mean);
	}
	return sqrt(temp / input.size());
}

double Graph::calculate_stdv_sample(double mean, vector<double>& input)//sample standard deviation
{
	double temp = 0;
	for (int i = 0; i < input.size(); i++)
	{
		temp += (input[i] - mean) * (input[i] - mean);
	}
	return sqrt(temp / (input.size() - 1));
}

Graph::Graph(int maxn, int maxl)
{
	maxnodes = maxn;
	maxlinks = maxl;
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

double Graph::dijkstra_independent(int source, int target)
{
	//Link** previous;
	//previous = new Link*[maxnodes];//de beste voorganger opslaan voor op het einde de optimale route uit te lezen
	BinaryMinHeap PQ(maxnodes);//aanmaken priority queue
	for (int i = 0; i < maxnodes; ++i)//score correct initialiseren
	{
		if (i == source)
		{
			n[i].score = 0;
			PQ.insert(&n[i]);
		}
		else
		{
			n[i].score = DBL_MAX;
			PQ.insert(&n[i]);
		}
	}
	Node* dn;
	dn = PQ.extractMin();
	while (dn->id != n[target].id)
	{
		double start_score = dn->score;//haalt kleinste score op
		//cout<<"tracking from"<<dn->id<<endl;
		Link* dl;
		dl = dn->first;//haalt een link address verbonden met de node op
		while (dl)
		{
			double score = start_score;
			//cout<<"startscore"<<start_score<<endl;
			score += (dl->optimaltt);//bereken de nieuwe score (reistijd in uur) van om naar die node te gaan
			if (score < n[dl->to].score)//vergelijk of die score kleiner is dan de huidige score die nu op de tweede node zit
			{
				//previous[dl->to]=dl;
				n[dl->to].score = score;//update de score
				//cout<<dl->to<<"/"<<n[dl->to].id<<"score"<<score<<endl;
				PQ.siftUp(n[dl->to].hdex);//zet de heap weer op orde je geeft de heapindex van de node mee
			}
			dl = dl->next;//haal een link address verbonden met de startnode op
		}//end while links
		dn = PQ.extractMin();
	}//end while not on target node
	return n[target].score;//is de traveltime
}

vector<double> Graph::dijkstra_independent_to_all_threaded(int source, vector<int>targets, int thread)
{
	BinaryMinHeaps PQp(maxnodes, thread);
	for (int i = 0; i < maxnodes; ++i)//score correct initialiseren
	{
		if (i == source)
		{
			np[i].score[thread] = 0;
			np[i].istarget[thread] = false;
			PQp.insert(&np[i], thread);
		}
		else
		{
			np[i].score[thread] = DBL_MAX;
			np[i].istarget[thread] = false;
			PQp.insert(&np[i], thread);
		}
	}
	for (int i = 0; i < (int)targets.size(); ++i)
	{
		np[targets[i]].istarget[thread] = true;
	}
	//targets might not be unique, so count how many unique nodes need to be found
	int targetgoal = 0;
	for (int i = 0; i < maxnodes; ++i)
	{
		if (np[i].istarget[thread] == true)
			++targetgoal;
	}
	int targetcount = 0;
	Nodep* dn;
	dn = PQp.extractMintwee(thread);
	while (targetcount < targetgoal)
	{
		if (dn->istarget[thread] == true)
			++targetcount;
		double start_score = dn->score[thread];//haalt kleinste score op
		//cout<<"tracking from"<<dn->id<<endl;
		Link* dl;
		dl = dn->first;//haalt een link address verbonden met de node op
		while (dl)
		{
			double score = start_score;
			score += (dl->optimaltt);//bereken de nieuwe score (reistijd in uur) van om naar die node te gaan
			if (score < np[dl->to].score[thread])//vergelijk of die score kleiner is dan de huidige score die nu op de tweede node zit
			{
				np[dl->to].score[thread] = score;//update de score
				//cout<<dl->to<<"/"<<np[dl->to].id<<"score"<<score<<endl;
				PQp.siftUp(np[dl->to].hdex[thread], thread);//zet de heap weer op orde je geeft de heapindex van de node mee
			}
			dl = dl->next;//haal een link address verbonden met de startnode op
		}//end while links
		dn = PQp.extractMintwee(thread);
	}//end while not on target node
	PQp.free(thread);//enkel de huidige thread clearen
	vector<double> targett;
	for (int t = 0; t < targets.size(); ++t)
	{
		targett.push_back(np[targets[t]].score[thread]);
	}
	return targett;
}

double Graph::dijkstra_dependent(int source, int target, double currenttime)//time-dependent dijkstra
{
	//Link** previous;
	//previous = new Link*[maxnodes];//de beste voorganger opslaan voor op het einde de optimale route uit te lezen
	BinaryMinHeap PQ(maxnodes);//aanmaken priority queue
	for (int i = 0; i < maxnodes; ++i)//score correct initialiseren
	{
		if (i == source)
		{
			n[i].score = currenttime;
			PQ.insert(&n[i]);
		}
		else
		{
			n[i].score = DBL_MAX;
			PQ.insert(&n[i]);
		}
	}
	Node* dn;
	dn = PQ.extractMin();
	while (dn->id != n[target].id)
	{
		double start_score = dn->score;//haalt kleinste score op
		//cout<<"tracking from"<<dn->id<<endl;
		Link* dl;
		dl = dn->first;//haalt een link address verbonden met de node op
		while (dl)
		{
			double score = start_score;
			//cout<<"startscore"<<start_score<<endl;
			int k = 0;
			while (score > time_periods[k + 1])
			{
				++k;
			}
			double interpol = dl->traveltime[k] + (double(dl->traveltime[k + 1] - dl->traveltime[k]) / (time_periods[k + 1] - time_periods[k])) * (score - time_periods[k]);
			score += interpol;
			if (score < n[dl->to].score)//vergelijk of die score kleiner is dan de huidige score die nu op de tweede node zit
			{
				//previous[dl->to]=dl;
				n[dl->to].score = score;//update de score
				//cout<<dl->to<<"/"<<n[dl->to].id<<"score"<<score<<endl;
				PQ.siftUp(n[dl->to].hdex);//zet de heap weer op orde je geeft de heapindex van de node mee
			}
			dl = dl->next;//haal een link address verbonden met de startnode op
		}//end while links
		dn = PQ.extractMin();
	}//end while not on target node
	return n[target].score - currenttime;//is de traveltime
}

vector<double> Graph::dijkstra_dependent_to_all_threaded(int source, vector<int>targets, double currenttime, int thread)
{
	BinaryMinHeaps PQp(maxnodes, thread);
	for (int i = 0; i < maxnodes; ++i)//score correct initialiseren
	{
		if (i == source)
		{
			np[i].score[thread] = currenttime;
			PQp.insert(&np[i], thread);
		}
		else
		{
			np[i].score[thread] = DBL_MAX;
			np[i].istarget[thread] = false;
			PQp.insert(&np[i], thread);
		}
	}
	for (int i = 0; i < (int)targets.size(); ++i)
	{
		np[targets[i]].istarget[thread] = true;
	}
	//targets might not be unique, so count how many unique nodes need to be found
	int targetgoal = 0;
	for (int i = 0; i < maxnodes; ++i)
	{
		if (np[i].istarget[thread] == true)
			++targetgoal;
	}
	int targetcount = 0;
	Nodep* dn;
	dn = PQp.extractMintwee(thread);
	while (targetcount < targetgoal)
	{
		if (dn->istarget[thread] == true)
			++targetcount;
		double start_score = dn->score[thread];//haalt kleinste score op
		//cout<<"tracking from"<<dn->id<<endl;
		Link* dl;
		dl = dn->first;//haalt een link address verbonden met de node op
		while (dl)
		{
			double score = start_score;
			unsigned int k = 0;
			while (score > time_periods[k + 1])
			{
				++k;
			}
			double interpol = dl->traveltime[k] + (double(dl->traveltime[k + 1] - dl->traveltime[k]) / (time_periods[k + 1] - time_periods[k])) * (score - time_periods[k]);
			score += interpol;
			if (score < np[dl->to].score[thread])//vergelijk of die score kleiner is dan de huidige score die nu op de tweede node zit
			{
				np[dl->to].score[thread] = score;//update de score
				PQp.siftUp(np[dl->to].hdex[thread], thread);//zet de heap weer op orde je geeft de heapindex van de node mee
			}
			dl = dl->next;//haal een link address verbonden met de startnode op
		}//end while links
		dn = PQp.extractMintwee(thread);
	}//end while not on target node
	PQp.free(thread);//enkel de huidige thread clearen
	vector<double> targett;
	for (int t = 0; t < targets.size(); ++t)
	{
		targett.push_back(np[targets[t]].score[thread] - currenttime);
	}
	return targett;
}


