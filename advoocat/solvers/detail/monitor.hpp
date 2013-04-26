#pragma once

#if defined(__linux__)
#  include <sys/prctl.h>

void monitor(float frac)
{
  char name[16];
  prctl(PR_GET_NAME, name, 0, 0, 0);
  // ...  1  0  0  % \0
  // ... 12 13 14 15 16
  int len = sprintf(&name[12], "%03d%%", int(frac * 100));
  prctl(PR_SET_NAME, name, 0, 0, 0);
}
#endif
