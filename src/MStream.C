#include <string.h>
#include "Communicate.h"
#include "MStream.h"

MIStream::MIStream(Communicate *c, int p, int t)
{
  cobj = c;
  PE = p;
  tag = t;
  msg = (StreamMessage *) 0;
  early = (StreamMessage *) 0;
  currentIndex = 0;
}

MIStream::~MIStream()
{
  if(msg!=0)
    CmiFree(msg);
}

MOStream::MOStream(Communicate *c, int p, int t, unsigned int size)
{
  cobj = c;
  PE = p;
  tag = t;
  bufLen = size;
  msgBuf = (StreamMessage *)CmiAlloc(sizeof(StreamMessage)+size);
  msgBuf->PE = CmiMyPe();
  msgBuf->tag = tag;
  msgBuf->len = 0;
  msgBuf->isLast = 0;
  msgBuf->index = 0;
  msgBuf->next = (StreamMessage *)0;
}

MOStream::~MOStream()
{
  if(msgBuf != 0)
    CmiFree(msgBuf);
}

MIStream *MIStream::Get(char *buf, int len)
{
  while(len) {
    if(msg==0) {
      if ( early && (early->index == currentIndex) ) {
        msg = early;
        early = early->next;
        msg->next = (StreamMessage *)0;
      } else {
        msg = (StreamMessage *) cobj->getMessage(PE, tag);
      }
      while ( msg->index != currentIndex ) {
        if ( (! early) || (early->index > msg->index) ) {
          msg->next = early;
          early = msg;
        } else {
          StreamMessage *cur = early;
          while ( cur->next && (cur->next->index < msg->index) ) {
            cur = cur->next;
          }
          msg->next = cur->next;
          cur->next = msg;
        }
        msg = (StreamMessage *) cobj->getMessage(PE, tag);
      }
      currentPos = 0;
      currentIndex += 1;
    }
    if(currentPos+len <= msg->len) {
      memcpy(buf, &(msg->data[currentPos]), len);
      currentPos += len;
      len = 0;
    } else {
      int b = msg->len-currentPos;
      memcpy(buf, &(msg->data[currentPos]), b);
      len -= b;
      buf += b;
      currentPos += b;
    }
    if(currentPos == msg->len) {
      CmiFree(msg);
      msg = 0;
    }
  }
  return this;
}

MOStream *MOStream::Put(char *buf, size_t len)
{
  while(len) {
    if(msgBuf->len + len <= bufLen) {
      memcpy(&(msgBuf->data[msgBuf->len]), buf, len);
      msgBuf->len += len;
      len = 0;
    } else {
      int b = bufLen - msgBuf->len;
      memcpy(&(msgBuf->data[msgBuf->len]), buf, b);
      msgBuf->len = bufLen;
      cobj->sendMessage(PE, (void *)msgBuf, bufLen+sizeof(StreamMessage));
      msgBuf->len = 0;
      msgBuf->isLast = 0;
      msgBuf->index += 1;
      len -= b;
      buf += b;
    }
  }
  return this;
}

void MOStream::end(void)
{
  msgBuf->isLast = 1;
  cobj->sendMessage(PE,(void*)msgBuf,msgBuf->len+sizeof(StreamMessage));
  msgBuf->len = 0;
  msgBuf->isLast = 0;
  msgBuf->index += 1;
}

