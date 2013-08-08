#include <iostream>
#include <string>

#include "gtest/gtest.h"

GTEST_API_ int main(int argc, char **argv) 
{

#if defined(_DEBUG) && defined(WIN32)

  _CrtSetReportMode( _CRT_ERROR, _CRTDBG_MODE_DEBUG );
  _CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
  //_CrtSetBreakAlloc(1268386);

#endif

  testing::InitGoogleTest(&argc, argv);
  int gtestRst = RUN_ALL_TESTS();

  return gtestRst;
}
