#ifndef OWNERBOX_H
#define OWNERBOX_H

template <class Owner, class Data>
class Box;

template <class Owner, class Data>
class OwnerBox {
friend class Box<Owner,Data>;
public:
  OwnerBox(Owner *o, void (Owner::*fn)() );
  ~OwnerBox();
      
  void open(Data* d);

  void close(); 
  Box<Owner,Data> *checkOut(void);

  void checkIn(Box<Owner,Data> * box); 
  int isOpen() { return (openCount); };

private:
  Owner *owner;
  void (Owner::*callback)();

  Data* data;

  int numberUsers, openCount, closeCount;
};

#endif
