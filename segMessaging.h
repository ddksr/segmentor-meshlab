#ifndef SEGMESSAGING_H
#define SEGMESSAGING_H

#include <string.h>

#include "libsegmentor/common.h"
#include <Qt>
#include <QMessageBox>

class segMessaging: public Messaging
{
  QWidget* win;
 public:
  segMessaging(QWidget*);
  
  void info(const char*);
  void error(const char*);
  void debug(const char*);
};

#endif
