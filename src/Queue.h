#ifndef QUEUE_H
#define QUEUE_H

// Quick and dirty queue for small numbers of items.

template<class T> class ShortQueue {
  private:

    class Elem
    {
      public:

      Elem(T &d) : next(NULL), data(d) {;}
      Elem *next;
      T data;
    };

    Elem *next;

  public:

    ShortQueue(void) : next(NULL) {;}

    int empty(void) { return ! next; }

    T pop(void)
    {
      T data = next->data;
      Elem *nn = next->next;
      delete next;
      next = nn;
      return data;
    }

    void append(T &data)
    {
      Elem *nn;
      for( nn = next ; nn->next ; nn = nn->next );
      nn->next = new Elem(data);
    }

};

#endif
