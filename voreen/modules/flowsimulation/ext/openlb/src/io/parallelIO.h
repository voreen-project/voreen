/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Jonas Latt
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#ifndef PARALLEL_IO_H
#define PARALLEL_IO_H

#include <streambuf>
#include <istream>
#include <ostream>
#include <fstream>
#include <iostream>

#include "communication/mpiManager.h"

namespace olb {

class ParBuf : public std::streambuf {
public:
  typedef enum {normal, serial} Modes;
  // There are two modes for proceeding of parallel i/o:
  // - normal: the default mode. all nodes write data, and all
  //           nodes receive data; this mode requests broadcasting.
  // - serial: only node 0 writes and receives data.
public:
  ParBuf(std::streambuf* _originalBuf);
  Modes   getMode() const;
  void    setMode(Modes _mode);
protected:
  int_type overflow (int_type c) override;
  std::streamsize xsputn(const char* s, std::streamsize num) override;

  int_type uflow() override;
  int_type underflow() override;
  std::streamsize xsgetn (char* s, std::streamsize num) override;
private:
  std::streambuf* originalBuf;
  Modes      mode;
};

class olb_ofstream : public std::ostream {
public:
  olb_ofstream();
  explicit olb_ofstream(const char * filename,
                        openmode mode = out | trunc );
  ~olb_ofstream() override;

  std::streambuf* rdbuf() const;
  bool is_open();
  void open(const char* filename, openmode mode = out | trunc);
  void close();
private:
  std::filebuf fbuf;
  ParBuf  mybuf;
};

class olb_ifstream : public std::istream {
public:
  olb_ifstream();
  explicit olb_ifstream(const char * filename,
                        openmode mode = in );
  ~olb_ifstream() override;

  std::streambuf* rdbuf() const;
  bool is_open();
  void open(const char* filename, openmode mode = in);
  void close();
private:
  std::filebuf fbuf;
  ParBuf  mybuf;
};

class olb_fstream : public std::iostream {
public:
  olb_fstream();
  explicit olb_fstream(const char * filename,
                       openmode mode = in | out );
  ~olb_fstream() override;

  std::streambuf* rdbuf() const;
  bool is_open();
  void open(const char* filename, openmode mode = in | out);
  void close();
private:
  std::filebuf fbuf;
  ParBuf  mybuf;
};

///////////////////////////////////////////////////////////////////
// Class ParBuf
///////////////////////////////////////////////////////////////////

inline ParBuf::ParBuf(std::streambuf* _originalBuf)
  : originalBuf(_originalBuf), mode(normal)
{ }

inline std::streambuf::int_type
ParBuf::overflow (std::streambuf::int_type c)
{
  int_type returnVal = c;
  if (c != EOF) {
#ifdef PARALLEL_MODE_MPI
    if (singleton::mpi().isMainProcessor()) {
#endif
      returnVal = originalBuf->sputc((char)c);
#ifdef PARALLEL_MODE_MPI
    }
    if (mode==normal) {
      singleton::mpi().bCast(&returnVal, 1);
    }
#endif
  }
  return returnVal;
}

inline ParBuf::Modes
ParBuf::getMode() const
{
  return mode;
}

inline void
ParBuf::setMode(ParBuf::Modes _mode)
{
  mode = _mode;
}

inline std::streamsize
ParBuf::xsputn(const char* s, std::streamsize num)
{
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()) {
#endif
    return originalBuf->sputn(s,num);
#ifdef PARALLEL_MODE_MPI
  }
  else {
    return num;
  }
#endif
}

inline std::streambuf::int_type
ParBuf::uflow()
{
  int_type value;
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()) {
#endif
    value = originalBuf->sbumpc();
#ifdef PARALLEL_MODE_MPI
  }
  if (mode==normal) {
    singleton::mpi().bCast(&value, 1);
  }
#endif
  return value;
}

inline std::streambuf::int_type
ParBuf::underflow()
{
  int_type value;
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()) {
#endif
    value = originalBuf->sgetc();
#ifdef PARALLEL_MODE_MPI
  }
  if (mode==normal) {
    singleton::mpi().bCast(&value, 1);
  }
#endif
  return value;
}

inline std::streamsize
ParBuf::xsgetn (char* s, std::streamsize num)
{
  std::streamsize sizeRead=0;
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()) {
#endif
    sizeRead = originalBuf->sgetn(s, num);
#ifdef PARALLEL_MODE_MPI
  }
  if (mode==normal) {
    int intSizeRead = (int) sizeRead;
    singleton::mpi().bCast(&intSizeRead, 1);
    singleton::mpi().bCast(s, intSizeRead);
    sizeRead = (std::streamsize) intSizeRead;
  }
#endif
  return sizeRead;
}

///////////////////////////////////////////////////////////////////
// Class olb_ofstream
///////////////////////////////////////////////////////////////////

inline olb_ofstream::olb_ofstream() : std::ostream(nullptr), fbuf(), mybuf(&fbuf)
{
  this->init(&mybuf);
}

inline olb_ofstream::olb_ofstream(const char * filename, openmode mode)
  : std::ostream(nullptr), fbuf(), mybuf(&fbuf)
{
  init(&mybuf);
  open(filename, mode);
}

inline olb_ofstream::~olb_ofstream()
{ }

inline std::streambuf*
olb_ofstream::rdbuf() const
{
  return const_cast<ParBuf*>(&mybuf);
}

inline bool
olb_ofstream::is_open()
{
  return fbuf.is_open();
}

inline void
olb_ofstream::open(const char* filename, openmode mode)
{
  int ok;
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()) {
#endif
    ok = (bool) fbuf.open(filename, mode | ios_base::out);
#ifdef PARALLEL_MODE_MPI
  }
  singleton::mpi().bCast(&ok, 1);
#endif
  if (!ok) {
    this->setstate(ios_base::failbit);
  }
}

inline void
olb_ofstream::close()
{
  int ok;
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()) {
#endif
    ok = (bool) fbuf.close();
#ifdef PARALLEL_MODE_MPI
  }
  singleton::mpi().bCast(&ok, 1);
#endif
  if (!ok) {
    setstate(ios_base::failbit);
  }
}



///////////////////////////////////////////////////////////////////
// Class olb_ifstream
///////////////////////////////////////////////////////////////////

inline olb_ifstream::olb_ifstream() : std::istream(nullptr), fbuf(), mybuf(&fbuf)
{
  init(&mybuf);
}

inline olb_ifstream::olb_ifstream(const char * filename, openmode mode)
  : std::istream(nullptr), fbuf(), mybuf(&fbuf)
{
  init(&mybuf);
  open(filename, mode);
}

inline olb_ifstream::~olb_ifstream()
{ }

inline std::streambuf*
olb_ifstream::rdbuf() const
{
  return const_cast<ParBuf*>(&mybuf);
}

inline bool
olb_ifstream::is_open()
{
  return fbuf.is_open();
}


inline void
olb_ifstream::open(const char* filename, openmode mode)
{
  int ok;
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()) {
#endif
    ok = (bool) fbuf.open(filename, mode | ios_base::in);
#ifdef PARALLEL_MODE_MPI
  }
  singleton::mpi().bCast(&ok, 1);
#endif
  if (!ok) {
    this->setstate(ios_base::failbit);
  }
}

inline void
olb_ifstream::close()
{
  int ok;
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()) {
#endif
    ok = (bool) fbuf.close();
#ifdef PARALLEL_MODE_MPI
  }
  singleton::mpi().bCast(&ok, 1);
#endif
  if (!ok) {
    setstate(ios_base::failbit);
  }
}


///////////////////////////////////////////////////////////////////
// Class olb_fstream
///////////////////////////////////////////////////////////////////

inline olb_fstream::olb_fstream() : std::iostream(nullptr), fbuf(), mybuf(&fbuf)
{
  this->init(&mybuf);
}

inline olb_fstream::olb_fstream(const char * filename, openmode mode)
  : std::iostream(nullptr), fbuf(), mybuf(&fbuf)
{
  init(&mybuf);
  open(filename, mode);
}

inline olb_fstream::~olb_fstream()
{ }

inline std::streambuf*
olb_fstream::rdbuf() const
{
  return const_cast<ParBuf*>(&mybuf);
}

inline bool
olb_fstream::is_open()
{
  return fbuf.is_open();
}

inline void
olb_fstream::open(const char* filename, openmode mode)
{
  int ok;
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()) {
#endif
    ok = (bool) fbuf.open(filename, mode);
#ifdef PARALLEL_MODE_MPI
  }
  singleton::mpi().bCast(&ok, 1);
#endif
  if (!ok) {
    this->setstate(ios_base::failbit);
  }
}

inline void
olb_fstream::close()
{
  int ok;
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()) {
#endif
    ok = (bool) fbuf.close();
#ifdef PARALLEL_MODE_MPI
  }
  singleton::mpi().bCast(&ok, 1);
#endif
  if (!ok) {
    setstate(ios_base::failbit);
  }
}

} // namespace olb

#endif
