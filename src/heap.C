class manheap;
class maxHeap;
#include "heap.h"
#include "iostream.h"
// Heap of pointers. The field to be compared is:

minHeap::minHeap(int size)
{
  h = new heapRecord[size];
  count = 0;
}

minHeap::numElements()
{
  return count;
}

void minHeap::insert(InfoRecord *x)
{
  h[count].info = x;
  h[count].deleted = 0;

  int current = count;
  count++;

  int parent = (current - 1)/2;
  while (current != 0)
    {
      if (h[current].info->load < h[parent].info->load)
	{
	  swap(current, parent);
	  current = parent;
	  parent = (current-1)/2;
	}
      else
	break;
    }
}

InfoRecord *minHeap::deleteMin()
{
  if (count == 0) return 0;

  InfoRecord *tmp = h[0].info;
  int best;

  h[0] = h[count-1];
  count--;

  int current = 0;
  int c1 = 1;
  int c2 = 2;
  while (c1 < count)
    {
      if (c2 >= count)
	best = c1;
      else
	{
	  if (h[c1].info->load < h[c2].info->load)
	    best = c1;
	  else
	    best = c2;
	}
      if (h[best].info->load < h[current].info->load)
	{
	  swap(best, current);
	  current = best;
	  c1 = 2*current + 1;
	  c2 = c1 + 1;
	}
      else
	break;
    }
  return tmp;
}

InfoRecord *minHeap::iterator(heapIterator *iter){
  iter->next = 1;
  if (count == 0)
    return 0;
  return h[0].info;
}
InfoRecord *minHeap::next(heapIterator *iter){
  if (iter->next >= count)
    return 0;
  iter->next += 1;
  return h[iter->next - 1].info;
}

//*****************


maxHeap::maxHeap(int size)
{
  h = new heapRecord[size];
  count = 0;
}

maxHeap::numElements()
{
  return count;
}

void maxHeap::insert(InfoRecord *x)
{
  h[count].info = x;
  h[count].deleted  = 0;
  int current = count;
  count++;

  int parent = (current - 1)/2;
  while (current != 0)
    {
      if (h[current].info->load > h[parent].info->load)
	{
	  swap(current, parent);
	  current = parent;
	  parent = (current-1)/2;
	}
      else
	break;
    }
}

InfoRecord *maxHeap::deleteMax()
{
  if (count == 0) return 0;
  InfoRecord *tmp = h[0].info;
  int best;

  h[0] = h[count-1];
  count--;

  int current = 0;
  int c1 = 1;
  int c2 = 2;
  while (c1 < count)
    {
      if (c2 >= count)
	best = c1;
      else
	{
	  if (h[c1].info->load > h[c2].info->load)
	    best = c1;
	  else
	    best = c2;
	}
      if (h[best].info->load > h[current].info->load)
	{
	  swap(best, current);
	  current = best;
	  c1 = 2*current + 1;
	  c2 = c1 + 1;
	}
      else
	break;
    }
  return tmp;
}


InfoRecord *maxHeap::iterator(heapIterator *iter){
  iter->next = 1;
  if (count == 0)
    return 0;
  return h[0].info;
}

InfoRecord *maxHeap::next(heapIterator *iter){
  if (iter->next >= count)
    return 0;
  iter->next += 1;
  return h[iter->next - 1].info;
}
