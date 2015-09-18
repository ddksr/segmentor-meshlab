#ifndef SEGRECOVER_DIALOG_H
#define SEGRECOVER_DIALOG_H

#include <QDialog>
#include <QDockWidget>
#include <qfiledialog.h>

#include <common/interfaces.h>
#include <meshlab/glarea.h>

#include "libsegmentor/common.h"
#include "ui_edit_segmentor_recover.h"

class segRecoverDialog : public QDockWidget
{
  Q_OBJECT

 public:
  segRecoverDialog(QWidget *);


 public:
  Ui::SegmentorRecoverDialog ui;
  void obtainSettings();
  void storeSettings();
  void recoverSettings();

  void closeEvent ( QCloseEvent * event ) ;

 private:
  QSettings *iniConfig;
  RecoverySettings config;
  
 signals:
  void closing();

};
		

#endif
