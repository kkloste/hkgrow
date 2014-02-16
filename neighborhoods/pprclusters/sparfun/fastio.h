/** @file fastio.h
 * Faster single-threaded IO operations.  
 * New versions of Visual Studio 2008 implement multi-threaded readers
 * with synchronized access to file handles and iostream objects.  For
 * single large-read operations, these synchronization mechanisms are not
 * a problem, but for many small-read operations, they incur significant 
 * overhead.  For example:
 *
 * The following codes:
 * 1) int s=0,k; while (n-->0) { fread(&k, sizeof(int), 1, f); s+=k }
 * 2) int *buf[n]; fread(&buf, sizeof(int), n, f); while(n-->0) { s+=*buf++; }
 * compute (or are intended to) compute identical values, but 1) takes 100x 
 * as long as 2) for n = 100000.
 *
 * To get the behavior of 1) with better performance, here are a few fastio
 * routines.
 */

/*
 * David F. Gleich
 * Copyright, Microsoft Corporation, 2008
 */

#ifndef FASTIO_H
#define FASTIO_H

#pragma warning(push)
#pragma warning(disable:4996)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <algorithm>
#include <inttypes.h>

#if defined(_WIN32) || defined(_WIN64)
#define ftell64 _ftelli64
#define fseek64 _fseeki64
#else
#define ftell64 ftell
#define fseek64 fseek
#endif

/** Implement a buffered reader for a single thread.  This class
 * avoids the thread-exclusion tools for reading files byte-by-byte
 * in recent C++/stdio implementations.
 */
class buffered_reader {
private:
    FILE *f;
    const size_t bufsize;
    size_t bufpos, bufused;
    unsigned char* buf;
    bool managed_file;

private:
    /** Refill the buffer from the current point in the file.
     */
    inline void refill() {
        bufused = fread(buf, sizeof(unsigned char), bufsize, f);
        bufpos = 0;
    }

    /** Refill the buffer so that we have bytes available in the buffer.
     * This call updates the buffer from it's current position so that
     * bytes of memory are stored in the buffer.
     *
     * It's used when 
     */
    inline void refill_tail(size_t bytes) {
        size_t tokeep = bufused - bufpos;
        memcpy(buf, &buf[bufpos], tokeep);
        bufused = tokeep;
        bufused += fread(&buf[tokeep],sizeof(unsigned char),bufsize-tokeep,f);
        bufpos = 0;
    }

public:
    buffered_reader(FILE *f_, size_t bufsize_)
        : f(f_), bufsize(bufsize_), bufpos(0), bufused(0), managed_file(false)
    {
        buf = (unsigned char*) malloc(bufsize);
    }

    buffered_reader(const char* filename, size_t bufsize_)
        : bufsize(bufsize_), bufpos(0), bufused(0), managed_file(true)
    {
        f = fopen(filename, "rbS");
        buf = (unsigned char*) malloc(bufsize);
    }

    buffered_reader(const buffered_reader& r) 
      : bufsize(r.bufsize), bufpos(r.bufpos), bufused(r.bufused),
        managed_file(r.managed_file), f(r.f)
    {
      assert(r.managed_file == false);
      buf = (unsigned char*) malloc(bufsize);
    }

    ~buffered_reader() { 
        free(buf);
        if (managed_file) fclose(f); 
    }

    int eof(void) {
        if (bufpos == bufused) {
            return feof(f); 
        } else {
            return 0;
        }
    }

    /** Get the size of this file
     *
     * This operation doesn't kill the buffer.
     */
    long long size() {
        int64_t pos = ftell64(f);
        if (pos < 0) {
            return (-1);
        }
        if (fseek(f, 0, SEEK_END) == 0) {
            int64_t size = ftell64(f);
            if (size >= 0 && fseek64(f, pos, SEEK_SET) == 0) {
                return (long long)size;
            } 
        }
        return (-1);
    }

    /** Seek to an arbitrary location in the file.
    * This operation clears the buffer.  
    * Use seekf for a buffer-preserving seek.
    */
    int seek(long long offset, int origin) {
        if (origin==SEEK_CUR) { offset -= (bufused-bufpos); }
        int rval=fseek64(f, offset, origin);
        bufpos = bufused; // this will cause the buffer to refill
        return (rval);
    }

    /** Advance to a new forward location in the file.
     * For small values of advance, this will not invalidate
     * the buffer and may be faster. If the advance position
     * is outside the buffer, it'll work just like seek.
     */
    int seekf(size_t advance) {
        if (bufpos + advance>=bufused) {
            return seek(advance, SEEK_CUR);
        } else {
            bufpos += advance; 
            return (0);
        }
    }

    /** Read data from the file.
     */
    size_t read(char *dst, size_t bytes) {
        size_t inbuf = 0, readbytes = 0;
        inbuf = std::min(bufused-bufpos,bytes);
        memcpy(dst, &buf[bufpos], inbuf);
        bytes -= inbuf;
        readbytes += inbuf;
        bufpos += inbuf;
        if (bytes > bufsize) {
            readbytes += fread(&dst[readbytes], sizeof(char), bytes, f);
            refill();
        } else if (bytes > 0) {
            refill();
            inbuf = std::min(bufused-bufpos,bytes);
            memcpy(&dst[readbytes], &buf[bufpos], inbuf);
            bytes -= inbuf;
            readbytes += inbuf;
            bufpos += inbuf;
        }
        return readbytes;
    }

    /** Read an item set
     * Works just like fread
     * TODO Fix this to work correctly when size*count exceeds size_t size
     */
    size_t read(void *buffer, size_t size, size_t count) {
      size_t nread = read((char*)buffer, size*count);
      return nread/size;
    }

    /** Get a pointer to our interal buffer with just enough bytes.
     *
     * This call is identical to read, but we give you access to 
     * the internal buffer instead of copying it.
     */
    bool bufptr(size_t bytes, const void** mem) {
        bool rval = largebufptr(bytes, mem);
        if (rval) { bufpos += bytes; }
        return rval;
    }

    /** Get a pointer to our interal buffer with just enough bytes.
     *
     * @return NULL on error, otherwise a pointer to the internal buffer.
     */
    const void* bufptr(size_t bytes) {
        const void* rval = NULL;
        if (bufptr(bytes, &rval)) { return (rval); }
        return (NULL);
    }

    /** Get a pointer to our interal buffer with too many bytes.
     * 
     * This call gives you a pointer with at least as many bytes as
     * requested but does not update the internal position.  Subsequent
     * operations use the previous position.  The idea is to get access
     * to more memory than you might need, but then only advance based
     * on how much you used.
     */
    bool largebufptr(size_t bytes, const void** mem) {
        if (bytes > bufsize) {
            // we cannot use this technique for very large reads
            return false;
        } else {
            if (bufused - bufpos < bytes) {
                // we can copy everything in, but we
                // have to refill the buffer in a special way
                refill_tail(bytes);
            }

            if (bufused - bufpos >= bytes) {
                *mem = &buf[bufpos];
                return (true);
            }
        }
        // we still didn't fit anything, 
        // this likely indicates end of file
        return (false);
    }

    /** Get a pointer to our interal buffer with too many bytes.
     *
     * @return NULL on error, otherwise a pointer to the internal buffer.
     */
    const void* largebufptr(size_t bytes) {
        const void* rval = NULL;
        if (largebufptr(bytes, &rval)) { return (rval); }
        return (NULL);
    }

    /** Read a particular type. 
     */
    inline int read() {
        int rval = 0;
        rval |= byte();
        rval |= byte() << 8;
        rval |= byte() << 16;
        rval |= byte() << 24;
        return (rval);
    }

    inline int byte() {
        if (bufpos >= bufused) {
            refill();
        }
        return buf[bufpos++];
    }
};

#pragma warning(pop)

#endif /* FASTIO_H */
