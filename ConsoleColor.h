// ConsoleColor.h

#pragma once
#include <ostream>

namespace term {

	enum class Color : int {
		reset = 0,
		red = 31,
		green = 32,
		yellow = 33,
		blue = 34,
		white = 37,
	};

	struct Set {
		int code;
	};

	inline std::ostream& operator<<(std::ostream& os, Set s) {
		return os << "\x1b[" << s.code << 'm';
	}

	inline Set fg(Color c) { return Set{ static_cast<int>(c) }; }
	inline Set reset() { return Set{ 0 }; }

} // namespace term