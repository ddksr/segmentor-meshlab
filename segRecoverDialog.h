#ifndef SEGRECOVER_DIALOG_H
#define SEGRECOVER_DIALOG_H

#include <QDialog>
#include <QSettings>
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
  segRecoverDialog(QWidget *, QSettings *);


 public:
  Ui::SegmentorRecoverDialog ui;

  void closeEvent ( QCloseEvent * event ) ;

 private:
  QSettings *iniConfig;
  RecoverySettings *config;
  
  void obtainSettings();
  void storeSettings();
  void recoverSettings();
  
 signals:
  void closing();

};
		

#endif
