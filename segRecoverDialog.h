#ifndef SEGRECOVER_DIALOG_H
#define SEGRECOVER_DIALOG_H

#include <QDialog>
#include <QDockWidget>
#include <QProgressBar>
#include <qfiledialog.h>

#include <common/interfaces.h>
#include <meshlab/glarea.h>

#include "segMesh.h"
#include "libsegmentor/common.h"
#include "libsegmentor/segmentor.h"
#include "ui_edit_segmentor_recover.h"

class ProgressBar : public ProgressIndicator {
  QProgressBar *bar;
  QLabel *label;
 public:
  ProgressBar(QProgressBar*, QLabel*);
  void setProcessName(const char*);
  void set(int);
  void clear();
  void clear(int, int);
};

class segRecoverDialog : public QDockWidget
{
  Q_OBJECT

 public:
  segRecoverDialog(QWidget *, MeshModel*, GLArea*);
  ~segRecoverDialog();


 public:
  Ui::SegmentorRecoverDialog ui;
  void obtainSettings();
  void storeSettings();
  void recoverSettings();

  void closeEvent ( QCloseEvent * event ) ;

 private:
  ProgressBar *pb;
  Segmentor *seg;
  QString settingsFile;
  RecoverySettings *config;

 private slots:
  void handleStoreSettings();
  void handleRestoreSettings();

  void placeSeeds();
  void grow();
  void selection();
  void finalSelection();
  
 signals:
  void closing();
  void released();

};

#endif