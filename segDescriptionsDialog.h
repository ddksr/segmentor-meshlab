#ifndef SEGDESCRIPTIONS_DIALOG_H
#define SEGDESCRIPTIONS_DIALOG_H

#include <QDialog>
#include <QDockWidget>
#include <QProgressBar>

#include "libsegmentor/common.h"
#include "libsegmentor/segmentor.h"
#include "ui_edit_segmentor_descriptions.h"

class segDescriptionsDialog : public QDockWidget
{
  Q_OBJECT

 public:
  segDescriptionsDialog(QWidget *);
  //~segDescriptionsDialog();

  void closeEvent (QCloseEvent * );

  Ui::SegmentorDescriptionsDialog ui;

 private:
  Segmentor *seg;

  void fillList();

 signals:
  void closing();
  void released();

};

#endif
