#if !defined(ASSERT_HPP)
  #define ASSERT_HPP

  #if defined(_DEBUG)
    void my_assert(const char* Condition, const char* FileName, int LineNumber);
    #define ASSERT(E) if (!(E))  my_assert(#E, __FILE__, __LINE__);
    #define VERIFY(E) if (!(E))  my_assert(#E, __FILE__, __LINE__);
  #else
    #define ASSERT(E)
    #define VERIFY(E) E;
  #endif

#endif
