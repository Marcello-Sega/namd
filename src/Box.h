#ifndef BOX_H
#define BOX_H

#include "ckdefs.h"
#include "char.h"
#include "c++interface.h"

class BoxOwner
{
private:
	boxIsClosed(void *) = 0;

friend OwnerBox
};

template T
class OwnerBox<class T>
{
public:
	OwnerBox(BoxOwner *o) :
		owner(o),
		number_of_users(0),
		accesses_remaining(0),
		data(0)
		{};
	void open(T* d, void *box_id = 0)
	{
		accesses_remaining = number_of_users;
		data = d;
	}
	void register(void) { ++number_of_users; }
	void unregister(void)
	{
		if ( ! number_of_users-- )
		{
			cprintf("OwnerBox::unregister() - no registrants remaining");
			number_of_users = 0;
		}
	}
	void unregisterAll(void) { number_of_users = 0; }

private:
	BoxOwner *owner;
	void *id;
	int number_of_users;
	int accesses_remaining;
	T* data;
	void close(void)
	{
		owner->boxIsClosed(id);
	}

friend Box<T>
};

template T
class Box<class T>
{
public:
	Box(OwnerBox<T>* o) : owner(o) {};
	T* open(void) { return owner->data; }
	void close(void)
	{
		if ( ! --(owner->accesses_remaining) ) owner->close();
	}

private:
	OwnerBox<T> *owner;	
};

#endif // BOX_H
