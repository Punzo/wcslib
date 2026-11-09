/*============================================================================
  WCSLIB 8.10 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2026, Mark Calabretta

  This file is part of WCSLIB.

  WCSLIB is free software: you can redistribute it and/or modify it under the
  terms of the GNU Lesser General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option)
  any later version.

  WCSLIB is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
  more details.

  You should have received a copy of the GNU Lesser General Public License
  along with WCSLIB.  If not, see http://www.gnu.org/licenses.

  Original author: Thomas Robitaille, The Astropy Project.
  $Id: twcsprintf_pthread.c,v 8.10 2026/09/27 05:44:34 mcalabre Exp $
*=============================================================================
* twcsprintf tests the inherent thread-safety of wcsprintf() when output is
* directed to an internal string buffer using thread-local storage (TLS).
*
* Each thread loops calling
*   wcsprintf_set(NULL);
*   wcsprintf("<thread-tid> <padding-string>\n");
*   const char *out = wcsprintf_buf();
*   -- check that `out` starts with the expected "<thread-tid>" prefix --
*
* The buffer contents read should always match what the thread just wrote.
* Without thread-local storage, the output and buffer resets from multiple
* threads interleave on the shared wcsprintf_buff and we may read back an
* empty string, the output from another thread, or a corrupted mixture.
* Moreover, concurrent buffer growth races realloc(), which may crash inside
* libc.
*
* To trigger the race condition, compile wcsprintf.c with WCSLIB_TLS set
* blank, which disables TLS (modify makedefs manually).
*
* Expected output:
*      with TLS: "OK: 0 corrupted out of <N> reads"
*   without TLS: many "corrupted" lines plus (often) a segfault
*
* An alternate to TLS would be to protect wcsprintf() and friends with a
* mutex - see commented-out code below.  This is less efficient than TLS
* because the threads have to wait on each other.  Moreover it does not
* provide an opportunity to aggregate messages from a thread before output,
* and requires additional code to create and manage the mutex.
*---------------------------------------------------------------------------*/

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wcsprintf.h"

#define NTHREADS    8
#define NITERATIONS 5000
#define PADDING     "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
#define PAD_BYTES   64  /* keep in sync with PADDING above */

// pthread_mutex_t mutex;

static unsigned long corruptions[NTHREADS];
static void *worker(void *arg);

int main(void)
{
  if (wcsprintf_set(NULL) != 0) {
    // Allocation failed.
    return 1;
  }

  // pthread_mutex_init(&mutex, NULL);

  pthread_t threads[NTHREADS];
  for (long t = 0; t < NTHREADS; t++) {
    if (pthread_create(&threads[t], NULL, worker, (void *)t) != 0) {
      perror("pthread_create");
      return 1;
    }
  }

  for (int t = 0; t < NTHREADS; t++) {
    pthread_join(threads[t], NULL);
  }

  // pthread_mutex_destroy(&mutex);

  unsigned long total_corruptions = 0;
  for (int t = 0; t < NTHREADS; t++) {
    total_corruptions += corruptions[t];
  }

  unsigned long total_reads = (unsigned long)NTHREADS * NITERATIONS;
  if (total_corruptions == 0) {
    printf("OK: 0 corrupted out of %lu reads\n", total_reads);
    return 0;
  } else {
    printf("RACE: %lu corrupted out of %lu reads\n",
           total_corruptions, total_reads);
    for (int t = 0; t < NTHREADS; t++) {
      printf("  thread %d: %lu\n", t, corruptions[t]);
    }
    return 1;
  }
}

//----------------------------------------------------------------------------
// Exercise wcsprintf_set() and wcsprintf() in a thread.
//----------------------------------------------------------------------------

static void *worker(void *arg)

{
  long tid = (long)arg;
  char prefix[32];

  // Our output starts with: "<tid> ".
  snprintf(prefix, sizeof(prefix), "%ld ", tid);
  size_t prefix_len = strlen(prefix);

  for (int i = 0; i < NITERATIONS; i++) {
    // pthread_mutex_lock(&mutex);

    // It wouldn't normally make sense in practice to invoke wcsprintf_set()
    // repeatedly like this.  Instead, it would be invoked once (outside the
    // loop), accumulate messages from the thread, print them if desired, then
    // clean up.  Note that if written directly to a file, say stdout,
    // messages from all threads would be interleaved.
    if (wcsprintf_set(NULL) != 0) {
      // Allocation failed, give up cleanly.
      // pthread_mutex_unlock(&mutex);
      return NULL;
    }

    wcsprintf("%ld %s\n", tid, PADDING);

    const char *out = wcsprintf_buf();
    if (out == NULL || strncmp(out, prefix, prefix_len) != 0) {
      corruptions[tid]++;
    }

    // pthread_mutex_unlock(&mutex);
  }

  // Release memory allocated for the buffer, thereby deleting the messages.
  // pthread_mutex_lock(&mutex);
  wcsprintf_set(stdout);
  // pthread_mutex_unlock(&mutex);

  return NULL;
}
