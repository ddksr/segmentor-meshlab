#include "segRecoverDialog.h"
#include "edit_segmentor_recover.h"

#include <Qt>
#include <QMessageBox>
#include <QDialog>

using namespace std;

void segRecoverDialog::closeEvent ( QCloseEvent * /*event*/ )
{
  emit closing();
}

segRecoverDialog::segRecoverDialog(QWidget *parent, QSettings *settings) : QDockWidget(parent)
{
  segRecoverDialog::ui.setupUi(this);

  this->setFeatures(QDockWidget::AllDockWidgetFeatures);
  this->setAllowedAreas(Qt::LeftDockWidgetArea);
  this->setFloating(true);

  iniConfig = settings;
  config = getDefaultRecoverySettings();
}

void segRecoverDialog::obtainSettings() {}

void segRecoverDialog::storeSettings() {}

void segRecoverDialog::recoverSettings() {}

