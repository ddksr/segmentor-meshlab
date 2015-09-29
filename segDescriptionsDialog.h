#ifndef SEGDESCRIPTIONS_DIALOG_H
#define SEGDESCRIPTIONS_DIALOG_H

#include <QDialog>
#include <QDockWidget>
#include <QProgressBar>
#include <QStringList>
#include <QStringListModel>

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
  void fillList();

  Ui::SegmentorDescriptionsDialog ui;

 private:
  Segmentor *seg;
  QStringList *stringList;
  QStringListModel *listModel;

  
  void setIJ(int, int&, int&);
  void reprepareDrawer();

 private slots:
  void handleGrowStep();
  void handleGrowFull();
  void handleDelete();
  void handleSelect();
  

 signals:
  void closing();
  void released();

};

#endif
