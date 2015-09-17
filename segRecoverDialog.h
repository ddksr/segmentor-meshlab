#ifndef SEGRECOVER_DIALOG_H
#define SEGRECOVER_DIALOG_H

#include <QDialog>
#include <QDockWidget>
#include <qfiledialog.h>

#include "ui_edit_segmentor_recover.h"

class segRecoverDialog : public QDockWidget
{
  Q_OBJECT

 private:
  
 public:
  segRecoverDialog(QWidget *parent);


 public:
  Ui::SegmentorRecoverDialog ui;

  void closeEvent ( QCloseEvent * event ) ;

 signals:
  void closing();

};
		

#endif
