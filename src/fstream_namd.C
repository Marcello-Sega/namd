
#include "fstream_namd.h"
#include "common.h"

void ofstream_namd::open(const char *_fname) {
  if ( fd ) NAMD_bug("ofstream_namd::open() called when file is already open");
  fname = _fname;
  fd = NAMD_open_text(_fname);
}

ofstream_namd& ofstream_namd::flush() {
    const std::string text = str();
    NAMD_write(fd, text.c_str(), text.size(), fname.c_str());
    str("");
    return *this;
}

void ofstream_namd::close() {
  flush();
  NAMD_close(fd, fname.c_str());
  fd = 0;
}

