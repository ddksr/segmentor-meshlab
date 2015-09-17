#include "segRecoverDialog.h"

#include <Qt>
#include <QMessageBox>
#include <QDialog>

using namespace std;

void segRecoverDialog::closeEvent ( QCloseEvent * /*event*/ )
{
  emit closing();
}

segRecoverDialog::segRecoverDialog(QWidget *parent) : QDockWidget(parent)
{
  segRecoverDialog::ui.setupUi(this);

  this->setFeatures(QDockWidget::AllDockWidgetFeatures);
  this->setAllowedAreas(Qt::LeftDockWidgetArea);
  this->setFloating(true);
}
