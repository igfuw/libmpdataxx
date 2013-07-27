#pragma once

#if defined(__linux__)
#  include <sys/prctl.h>
#endif

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      void monitor(float frac)
      {
#if defined(__linux__)
	char name[16];
	prctl(PR_GET_NAME, name, 0, 0, 0);
	static int len = strlen(name);
	// ...     1  0  0  % \0
	// ... 10 11 12 13 14 15
	sprintf(
	  &name[std::min(len, 10)], // taking care of short thread-names
	  " %3d%%", 
	  int(frac * 100)
	);
	prctl(PR_SET_NAME, name, 0, 0, 0);
#endif
      }
    };
  };
};
