#if !defined(ASSERT_HPP)
  #define ASSERT_HPP

  #if defined(_DEBUG)
    void assert(char* Condition, char* FileName, int LineNumber);
    #define ASSERT(E) if (!(E))  assert(#E, __FILE__, __LINE__);
    #define VERIFY(E) if (!(E))  assert(#E, __FILE__, __LINE__);
  #else
    #define ASSERT(E)
    #define VERIFY(E) E;
  #endif

#endif
