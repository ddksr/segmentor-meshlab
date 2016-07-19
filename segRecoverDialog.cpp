#include <Qt>
#include <QMessageBox>
#include <QSettings>
#include <QDir>
#include <QDialog>

#include "libsegmentor/common.h"
#include "libsegmentor/segmentor.h"
#include "segRecoverDialog.h"
#include "segDescriptionsDialog.h"
#include "edit_segmentor_recover.h"
#include "segMesh.h"
#include "segDrawer.h"

using namespace std;

void segRecoverDialog::closeEvent ( QCloseEvent * /*event*/ )
{
  emit closing();
}

segRecoverDialog::~segRecoverDialog() {
  delete config;
  delete mesh;
  delete d;
  delete pb;
  delete message;

  if (desc != NULL) {
	delete desc;
  }
}

segRecoverDialog::segRecoverDialog(QWidget *_parent, MeshModel* m, GLArea* gl) : QDockWidget(_parent)
{
  gla = gl;
  segRecoverDialog::ui.setupUi(this);

  this->setFeatures(QDockWidget::AllDockWidgetFeatures);
  this->setAllowedAreas(Qt::LeftDockWidgetArea);
  this->setFloating(true);

  QObject::connect(ui.btnSaveSettings, SIGNAL(released()), this, SLOT(handleStoreSettings()));
  QObject::connect(ui.btnRestoreSettings, SIGNAL(released()), this, SLOT(handleRestoreSettings()));
  

  settingsFile = QDir::homePath() + QDir::separator() + "segmentor.ini";

  config = new RecoverySettings{
	{false, 3.0, 2.0}, //plane
	{true, 6.0, 5.0}, // sq
	{false, 3.0, 2.0}, // surface2
	{false, 3.0, 2.0}, // sphere
	{false, 3.0, 2.0}, // cylinder
	{false, 3.0, 2.0}, // cone
	{false, 3.0, 2.0}, // torus
	
	{false, 3.0, 2.0}, // asq
	{false, 3.0, 2.0}, // tsq
	{false, 3.0, 2.0}, // bsq
	
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
	false,
	false
  };

  segRecoverDialog::recoverSettings();

  mesh = new segMesh(m);
  d = new MeshlabDrawer((image*) mesh);
  pb = new ProgressBar(ui.progressBar, ui.progressLabel);
  seg = Segmentor::Instance();
  message = new segMessaging(gl->window()->parentWidget());

  seg->setUp(config, mesh, (Drawer*)d, (ProgressIndicator*) pb, message);

  QObject::connect(ui.btnPlaceSeeds, SIGNAL(clicked()), this, SLOT(placeSeeds()));
  QObject::connect(ui.btnGrow, SIGNAL(clicked()), this, SLOT(grow()));
  QObject::connect(ui.btnSelection, SIGNAL(clicked()), this, SLOT(selection()));
  QObject::connect(ui.btnFinalSelection, SIGNAL(clicked()), this, SLOT(finalSelection()));
  QObject::connect(ui.btnShowDescriptions, SIGNAL(clicked()), this, SLOT(showDescriptions()));
  QObject::connect(ui.btnDeleteWrong, SIGNAL(clicked()), this, SLOT(deleteWrong()));

  QObject::connect(ui.btnMerging, SIGNAL(clicked()), this, SLOT(merging()));
  QObject::connect(ui.btnFillGaps, SIGNAL(clicked()), this, SLOT(fillGaps()));
  QObject::connect(ui.btnRefinement, SIGNAL(clicked()), this, SLOT(refinement()));
  QObject::connect(ui.btnBasicRas, SIGNAL(clicked()), this, SLOT(basicRAS()));
  QObject::connect(ui.btnImprovedRas, SIGNAL(clicked()), this, SLOT(improvedRAS()));
  QObject::connect(ui.btnClear, SIGNAL(clicked()), this, SLOT(clear()));

  desc = NULL;
  seg->enableMessaging(true);
}

void segRecoverDialog::obtainSettings() {
  config->k2 = ui.inputK2->text().toFloat();
  config->k3 = ui.inputK3->text().toFloat();
  config->planarNormalCompatibility = ui.inputPlanarNormalCompatibility->text().toFloat();
  config->growingSteps = ui.inputGrowingSteps->text().toInt();
  config->selections = ui.inputSelections->text().toInt();
  config->seedSize = ui.inputSeedSize->text().toInt();
  config->regionCompatibility = ui.inputRegionCompatibility->text().toInt();
  config->maxNumOfSeeds = ui.inputMaxNumOfSeeds->text().toInt();
  config->varianceOfNoise = ui.inputVarianceOfNoise->text().toFloat();
  config->postProcessingMaxError = ui.inputPostProcessingMaxError->text().toFloat();
  
  config->isNewDiscrepancy = ui.cbNewDiscrepancy->isChecked();
  config->usePostProcessing = ui.cbPostProcessing->isChecked();
  config->useStatistics = ui.cbUseStatistics->isChecked();
  config->useLookup = ui.cbUseLookup->isChecked();
  config->useAsqKf = ui.cbUseAsqKf->isChecked();

  config->plane.on = ui.cbPlane->isChecked();
  config->plane.dist = ui.inputPlaneDist->text().toFloat();
  config->plane.err = ui.inputPlaneErr->text().toFloat();

  config->superquadric.on = ui.cbSq->isChecked();
  config->superquadric.dist = ui.inputSqDist->text().toFloat();
  config->superquadric.err = ui.inputSqErr->text().toFloat();

  config->secondOrderSurface.on = ui.cb2os->isChecked();
  config->secondOrderSurface.dist = ui.input2osDist->text().toFloat();
  config->secondOrderSurface.err = ui.input2osErr->text().toFloat();

  config->sphere.on = ui.cbSphere->isChecked();
  config->sphere.dist = ui.inputSphereDist->text().toFloat();
  config->sphere.err = ui.inputSphereErr->text().toFloat();

  config->cylinder.on = ui.cbCylinder->isChecked();
  config->cylinder.dist = ui.inputCylinderDist->text().toFloat();
  config->cylinder.err = ui.inputCylinderErr->text().toFloat();

  config->cone.on = ui.cbCone->isChecked();
  config->cone.dist = ui.inputConeDist->text().toFloat();
  config->cone.err = ui.inputConeErr->text().toFloat();

  config->torus.on = ui.cbTorus->isChecked();
  config->torus.dist = ui.inputTorusDist->text().toFloat();
  config->torus.err = ui.inputTorusErr->text().toFloat();

  config->asq.on = ui.cbASq->isChecked();
  config->asq.dist = ui.inputASqDist->text().toFloat();
  config->asq.err = ui.inputASqErr->text().toFloat();
  
  config->tsq.on = ui.cbTSq->isChecked();
  config->tsq.dist = ui.inputTSqDist->text().toFloat();
  config->tsq.err = ui.inputTSqErr->text().toFloat();

  config->bsq.on = ui.cbBSq->isChecked();
  config->bsq.dist = ui.inputBSqDist->text().toFloat();
  config->bsq.err = ui.inputBSqErr->text().toFloat();

}

void segRecoverDialog::storeSettings() {
  segRecoverDialog::obtainSettings();
  QSettings iniConfig(settingsFile, QSettings::NativeFormat);
  
  iniConfig.setValue("recover/k2", (double)config->k2);
  iniConfig.setValue("recover/k3", (double)config->k3);
  iniConfig.setValue("recover/planarNormalCompatibility", (double)config->planarNormalCompatibility);
  iniConfig.setValue("recover/growingSteps", config->growingSteps);
  iniConfig.setValue("recover/selections", config->selections);
  iniConfig.setValue("recover/seedSize", config->seedSize);
  iniConfig.setValue("recover/regionCompatibility", (double)config->regionCompatibility);
  iniConfig.setValue("recover/maxNumOfSeeds", config->maxNumOfSeeds);
  iniConfig.setValue("recover/varianceOfNoise", (double)config->varianceOfNoise);
  iniConfig.setValue("recover/postProcessingMaxError", (double)config->postProcessingMaxError);

  iniConfig.setValue("recover/isNewDiscrepancy", config->isNewDiscrepancy);
  iniConfig.setValue("recover/usePostProcessing", config->usePostProcessing);
  iniConfig.setValue("recover/useStatistics", config->useStatistics);
  iniConfig.setValue("recover/useLookup", config->useLookup);
  iniConfig.setValue("recover/useAsqKf", config->useAsqKf);

  iniConfig.setValue("plane/on", config->plane.on);
  iniConfig.setValue("plane/maxDistance", (double)config->plane.dist);
  iniConfig.setValue("plane/maxError", (double)config->plane.err);

  iniConfig.setValue("superquadric/on", config->superquadric.on);
  iniConfig.setValue("superquadric/maxDistance",(double)config->superquadric.dist);
  iniConfig.setValue("superquadric/maxError", (double)config->superquadric.err);

  iniConfig.setValue("secondOrderSurface/on", config->secondOrderSurface.on);
  iniConfig.setValue("secondOrderSurface/maxDistance", (double)config->secondOrderSurface.dist);
  iniConfig.setValue("secondOrderSurface/maxError", (double)config->secondOrderSurface.err);

  iniConfig.setValue("sphere/on", config->sphere.on);
  iniConfig.setValue("sphere/maxDistance", (double)config->sphere.dist);
  iniConfig.setValue("sphere/maxError", (double)config->sphere.err);

  iniConfig.setValue("cylinder/on", config->cylinder.on);
  iniConfig.setValue("cylinder/maxDistance", (double)config->cylinder.dist);
  iniConfig.setValue("cylinder/maxError", (double)config->cylinder.err);

  iniConfig.setValue("cone/on", config->cone.on);
  iniConfig.setValue("cone/maxDistance", (double)config->cone.dist);
  iniConfig.setValue("cone/maxError", (double)config->cone.err);

  iniConfig.setValue("torus/on", config->torus.on);
  iniConfig.setValue("torus/maxDistance", (double)config->torus.dist);
  iniConfig.setValue("torus/maxError", (double)config->torus.err);

  iniConfig.setValue("asq/on", config->asq.on);
  iniConfig.setValue("asq/maxDistance", (double)config->asq.dist);
  iniConfig.setValue("asq/maxError", (double)config->asq.err);

  iniConfig.setValue("tsq/on", config->tsq.on);
  iniConfig.setValue("tsq/maxDistance", (double)config->tsq.dist);
  iniConfig.setValue("tsq/maxError", (double)config->tsq.err);

  iniConfig.setValue("bsq/on", config->bsq.on);
  iniConfig.setValue("bsq/maxDistance", (double)config->bsq.dist);
  iniConfig.setValue("bsq/maxError", (double)config->bsq.err);

}

void segRecoverDialog::recoverSettings() {
  QSettings iniConfig(settingsFile, QSettings::NativeFormat);
  
  ui.inputK2->setText(iniConfig.value("recover/k2", config->k2).toString());
  ui.inputK3->setText(iniConfig.value("recover/k3", config->k3).toString());
  ui.inputPlanarNormalCompatibility->setText(iniConfig.value("recover/planarNormalCompatibility", config->planarNormalCompatibility).toString());
  ui.inputGrowingSteps->setText(iniConfig.value("recover/growingSteps", config->growingSteps).toString());
  ui.inputSelections->setText(iniConfig.value("recover/selections", config->selections).toString());
  ui.inputSeedSize->setText(iniConfig.value("recover/seedSize", config->seedSize).toString());
  ui.inputRegionCompatibility->setText(iniConfig.value("recover/regionCompatibility", config->regionCompatibility).toString());
  ui.inputMaxNumOfSeeds->setText(iniConfig.value("recover/maxNumOfSeeds", config->maxNumOfSeeds).toString());
  ui.inputVarianceOfNoise->setText(iniConfig.value("recover/varianceOfNoise", config->varianceOfNoise).toString());
  ui.inputPostProcessingMaxError->setText(iniConfig.value("recover/postProcessingMaxError", config->postProcessingMaxError).toString());

  ui.cbNewDiscrepancy->setChecked(iniConfig.value("recover/isNewDiscrepancy", config->isNewDiscrepancy).toBool());
  ui.cbPostProcessing->setChecked(iniConfig.value("recover/usePostProcessing", config->usePostProcessing).toBool());
  ui.cbUseStatistics->setChecked(iniConfig.value("recover/useStatistics", config->useStatistics).toBool());
  ui.cbUseLookup->setChecked(iniConfig.value("recover/useLookup", config->useLookup).toBool());
  ui.cbUseAsqKf->setChecked(iniConfig.value("recover/useAsqKf", config->useAsqKf).toBool());

  ui.cbPlane->setChecked(iniConfig.value("plane/on", config->plane.on).toBool());
  ui.inputPlaneDist->setText(iniConfig.value("plane/maxDistance", config->plane.dist).toString());
  ui.inputPlaneErr->setText(iniConfig.value("plane/maxError", config->plane.err).toString());

  ui.cbSq->setChecked(iniConfig.value("superquadric/on", config->superquadric.on).toBool());
  ui.inputSqDist->setText(iniConfig.value("superquadric/maxDistance", config->superquadric.dist).toString());
  ui.inputSqErr->setText(iniConfig.value("superquadric/maxError", config->superquadric.err).toString());

  ui.cb2os->setChecked(iniConfig.value("secondOrderSurface/on", config->secondOrderSurface.on).toBool());
  ui.input2osDist->setText(iniConfig.value("secondOrderSurface/maxDistance", config->secondOrderSurface.dist).toString());
  ui.input2osErr->setText(iniConfig.value("secondOrderSurface/maxError", config->secondOrderSurface.err).toString());

  ui.cbSphere->setChecked(iniConfig.value("sphere/on", config->sphere.on).toBool());
  ui.inputSphereDist->setText(iniConfig.value("sphere/maxDistance", config->sphere.dist).toString());
  ui.inputSphereErr->setText(iniConfig.value("sphere/maxError", config->sphere.err).toString());

  ui.cbCylinder->setChecked(iniConfig.value("cylinder/on", config->cylinder.on).toBool());
  ui.inputCylinderDist->setText(iniConfig.value("cylinder/maxDistance", config->cylinder.dist).toString());
  ui.inputCylinderErr->setText(iniConfig.value("cylinder/maxError", config->cylinder.err).toString());

  ui.cbCone->setChecked(iniConfig.value("cone/on", config->cone.on).toBool());
  ui.inputConeDist->setText(iniConfig.value("cone/maxDistance", config->cone.dist).toString());
  ui.inputConeErr->setText(iniConfig.value("cone/maxError", config->cone.err).toString());

  ui.cbTorus->setChecked(iniConfig.value("torus/on", config->torus.on).toBool());
  ui.inputTorusDist->setText(iniConfig.value("torus/maxDistance", config->torus.dist).toString());
  ui.inputTorusErr->setText(iniConfig.value("torus/maxError", config->torus.err).toString());

  ui.cbASq->setChecked(iniConfig.value("asq/on", config->asq.on).toBool());
  ui.inputASqDist->setText(iniConfig.value("asq/maxDistance", config->asq.dist).toString());
  ui.inputASqErr->setText(iniConfig.value("asq/maxError", config->asq.err).toString());

  ui.cbTSq->setChecked(iniConfig.value("tsq/on", config->tsq.on).toBool());
  ui.inputTSqDist->setText(iniConfig.value("tsq/maxDistance", config->tsq.dist).toString());
  ui.inputTSqErr->setText(iniConfig.value("tsq/maxError", config->tsq.err).toString());

  ui.cbBSq->setChecked(iniConfig.value("bsq/on", config->bsq.on).toBool());
  ui.inputBSqDist->setText(iniConfig.value("bsq/maxDistance", config->bsq.dist).toString());
  ui.inputBSqErr->setText(iniConfig.value("bsq/maxError", config->bsq.err).toString());
}

void segRecoverDialog::handleStoreSettings() {
  segRecoverDialog::storeSettings();
  QMessageBox::information(this->parentWidget(), "Segmentor recovery", "Settings stored");
}

void segRecoverDialog::handleRestoreSettings() {
  segRecoverDialog::recoverSettings();
  QMessageBox::information(this->parentWidget(), "Segmentor recovery", "Settings restored");
}

void segRecoverDialog::placeSeeds() {
  qDebug() << "SLOT: place seeds";
  obtainSettings();

  if (ui.radioSelectedPtSeeds->isChecked()) {
	mesh->setSelectedPoints();
  } else if (mesh->numOfSelectedPoints > 0) {
	mesh->numOfSelectedPoints = 0;
	delete mesh->selectedPoints;
  }
  
  seg->placeSeeds();
  gla->update();
}

void segRecoverDialog::grow() {
  qDebug() << "SLOT: grow";
  obtainSettings();
  seg->grow();
  gla->update();
}

void segRecoverDialog::selection() {
  qDebug() << "SLOT: selection";
  obtainSettings();
  seg->selection();
  gla->update();
}

void segRecoverDialog::finalSelection() {
  qDebug() << "SLOT: finalSelection";
  obtainSettings();
  seg->finalSelection();
  gla->update();
}

void segRecoverDialog::deleteWrong() {
  obtainSettings();
  seg->deleteWrong();
  gla->update();
}


void segRecoverDialog::merging() {
  obtainSettings();
  seg->mergeDescriptions();
  gla->update();
}

void segRecoverDialog::fillGaps() {
  obtainSettings();
  seg->fillGaps();
  gla->update();
}
void segRecoverDialog::refinement() {
  obtainSettings();
  seg->intersectionRefinement();
  gla->update();
}

void segRecoverDialog::basicRAS() {
  seg->enableMessaging(false);
  obtainSettings();
  if (ui.radioSelectedPtSeeds->isChecked()) {
	mesh->setSelectedPoints();
  } else if (mesh->numOfSelectedPoints > 0) {
	mesh->numOfSelectedPoints = 0;
	delete mesh->selectedPoints;
  }
  seg->placeSeeds();
  gla->update();
  int nk;
  for (int i = 0; i < config->selections; i++) {
	nk = 1 << i;
	for (int j = 0; j < nk; j++) {
	  seg->grow();
	  gla->update();
	}
	seg->selection();
	gla->update();
  }
  message->info("Basic RAS completed");
  seg->enableMessaging(true);
}
void segRecoverDialog::improvedRAS() {
  seg->enableMessaging(false);
  obtainSettings();
  if (ui.radioSelectedPtSeeds->isChecked()) {
	mesh->setSelectedPoints();
  } else if (mesh->numOfSelectedPoints > 0) {
	mesh->numOfSelectedPoints = 0;
	delete mesh->selectedPoints;
  }
  seg->placeSeeds();

  int nk;
  for (int i = 0; i < config->selections; i++) {
	nk = 1 << i;
	for (int j = 0; j < nk; j++) {
	  seg->grow();
	  gla->update();
	}
	seg->selection();
  }
  seg->finalSelection();
  seg->fillGaps();
  int ok = 1;
  while(ok) {
	ok = 0;
	for (int i = 0; i < seg->numOfDescriptions; i++) { ok+= seg->descriptions[i].n; }
	seg->intersectionRefinement();
	seg->mergeDescriptions();
	for (int i = 0; i < seg->numOfDescriptions; i++) {
	  seg->descriptions[i].throw_away();
	  ok -= seg->descriptions[i].n;
	}
  }
  seg->enableMessaging(true);
  message->info("Improved RAS completed");
}

void segRecoverDialog::showDescriptions() {
  if (desc == NULL) {
	desc = new segDescriptionsDialog(this->parentWidget(), gla);
	qDebug() << "descriptions dialog created";
  } else {
	desc->fillList();
	qDebug() << "descriptions dialog exists";
  }
  desc->show();
}

void segRecoverDialog::clear() {
  seg->setUp(config, mesh, (Drawer*)d, (ProgressIndicator*) pb, message);
  seg->enableMessaging(true);
}

ProgressBar::ProgressBar(QProgressBar *_bar, QLabel *_label) {
  bar = _bar;
  label = _label;
  
  clear(0, 100);
}

void ProgressBar::set(int x) {
  bar->setValue(x);
}

void ProgressBar::clear() {
  clear(0, 100);
  label->setText(QString(""));
}

void ProgressBar::clear(int min, int max) {
  bar->reset();
  bar->setRange(min, max);
  label->setText(QString(""));
}

void ProgressBar::setProcessName(const char* name) {
  label->setText(QString(name));
}

void ProgressBar::inc() {
  bar->setValue(bar->value() + 1);
}
