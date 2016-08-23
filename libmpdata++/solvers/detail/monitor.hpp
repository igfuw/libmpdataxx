// inspired by pulseaudio/src/pulsecore/thread-posix.c

// TODO: the HAVE_PTHREAD_SETNAME_NP and HAVE_PTHREAD_GETNAME_NP are not yet defined anywhere!

#pragma once

#if defined(__linux__)
#  include <sys/prctl.h> // Linux ''standard''
#elif defined(HAVE_PTHREAD_SETNAME_NP) && defined(HAVE_PTHREAD_GETNAME_NP)
#  include <pthread.h>   // POSIX ''non-standard''
#endif

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      void monitor(float frac)
      {
	char name[17];
#if defined(__linux__)
        name[16] = '\0';
	prctl(PR_GET_NAME, name, 0, 0, 0);
#elif defined(HAVE_PTHREAD_SETNAME_NP) && defined(HAVE_PTHREAD_GETNAME_NP)
        pthread_getname_np(pthread_self(), name, 16);
#endif
	static int len = strlen(name);
	// ...     1  0  0  % \0
	// ... 10 11 12 13 14 15
	sprintf(
	  &name[std::min(len, 10)], // taking care of short thread-names
	  " %3d%%", 
	  int(frac * 100)
	);
#if defined(__linux__)
	prctl(PR_SET_NAME, name, 0, 0, 0);
#elif defined(HAVE_PTHREAD_SETNAME_NP) && defined(HAVE_PTHREAD_GETNAME_NP)
        pthread_setname_np(name);
#endif
      }
    }
  }
}
