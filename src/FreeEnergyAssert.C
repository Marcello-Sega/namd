#include <iostream.h>
#include "InfoStream.h"
#include "FreeEnergyAssert.h"

void my_assert(const char* Condition, const char* FileName, int LineNumber) {
  iout << endl << endi;
  iout << "Assertion: " << "(" << Condition << ")," << " failed" << endl << endi;
  iout << "   in: " << FileName << ", " << "line: " << LineNumber << endl << endi;
  iout << endl << endi;
}
