#include "segRecoverDialog.h"
#include "edit_segmentor_recover.h"

#include <Qt>
#include <QMessageBox>
#include <QSettings>
#include <QDir>

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

  QString settingsFile = QDir::homePath() + QDir::separator() + "segmentor.ini";
	
  QSettings settings(settingsFile, QSettings::NativeFormat);
  iniConfig = &settings;

  config = {
	{false, 3.0, 2.0},
	{false, 6.0, 5.0},
	{false, 3.0, 2.0},
	{false, 3.0, 2.0},
	{false, 3.0, 2.0},
	{false, 3.0, 2.0},
	{false, 3.0, 2.0},
	0.2,
	0.3,
	0.75, 
	8,
	5,
	16,
	0.95,
	100,
	1.0,
	0.1,
	true,
	false
  };

  segRecoverDialog::recoverSettings();
}

void segRecoverDialog::obtainSettings() {
  config.k2 = ui.inputK2->text().toFloat();
}

void segRecoverDialog::storeSettings() {
  segRecoverDialog::obtainSettings();
  
  iniConfig->setValue("recover/k2", config.k2);
  iniConfig->setValue("recover/k3", config.k3);
}

void segRecoverDialog::recoverSettings() {
  ui.inputK2->setText(iniConfig->value("recover/k2", config.k2).toString());
  ui.inputK3->setText(iniConfig->value("recover/k3", config.k3).toString());
}

