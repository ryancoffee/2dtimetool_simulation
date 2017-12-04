#include <iostream>

namespace DebugOps
{
	int return_early(const int r){
		std::cerr << "\n"
			<< "===========================================\n"
			<< "======== returning early ==================\n"
			<< "======== exit code " << r << "================\n"
			<< std::flush;
		return r;
	}
}
