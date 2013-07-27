#pragma once

// TODO: namespace detail

#if defined(__linux__)
#  include <sys/prctl.h>

void monitor(float frac)
{
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
}
#endif
