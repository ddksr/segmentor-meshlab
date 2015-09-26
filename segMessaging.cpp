#include <QDebug>
#include <QMessageBox>
#include "segMessaging.h"

segMessaging::segMessaging(QWidget* _win) {
  win = _win;
}

void segMessaging::info(const char* text) {
  QMessageBox::information(win, QString("Segmentor"), QString(text));
}

void segMessaging::error(const char* text) {
  QMessageBox::critical(win, QString("Segmentor"), QString(text));
}

void segMessaging::debug(const char* text) {
  qDebug() << text;
}
