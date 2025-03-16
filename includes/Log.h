#pragma once

#include <iostream>
#include <sstream>

class Log {
	public:
		void append_indices(std::ostringstream& oss) { }

		void append_indices(std::ostringstream& oss, int i) { oss << "[" << i << "]"; }

		template <typename... Args>
		void append_indices(std::ostringstream& oss, int i, Args... rest) {
			oss << "[" << i << "]";
			append_indices(oss, rest...);
		}

		template <typename T, typename... Args>
		void log_value(const std::string& var_name, const T& value, Args... indices) {
			std::ostringstream oss;
			oss << var_name;
			append_indices(oss, indices...);
			oss << " = " << value;
			std::cout << oss.str() << std::endl;
		}
};


