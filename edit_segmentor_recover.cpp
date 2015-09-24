#include <Qt>
#include <QSettings>

#include "edit_segmentor_recover.h"
#include "segRecoverDialog.h"
#include "segDrawer.h"

using namespace std;
using namespace vcg;

EditSegmentorRecoverPlugin::EditSegmentorRecoverPlugin() {
  qFont.setFamily("Helvetica");
  qFont.setPixelSize(12);
	
  qDebug() << "Init Edit Segmentor Recover";

  recoverDialog = NULL;

  qDebug() << "End Init Edit Segmentor Recover";
}

const QString EditSegmentorRecoverPlugin::Info() 
{
  return tr("Segmentor for Meshlab.");
}

bool EditSegmentorRecoverPlugin::StartEdit(MeshModel &_m, GLArea *_gla ) {
  this->model = &_m;
  gla = _gla;

  qDebug() << "SEG: Start Edit";

  if (recoverDialog != NULL) {
	//delete recoverDialog;
  }

  recoverDialog = new segRecoverDialog(gla->window(), model, gla);
  recoverDialog->show();
  
  editDialogOn = true;
  
  return true;
}

void EditSegmentorRecoverPlugin::EndEdit(MeshModel &_m, GLArea *_gla) {
  qDebug() << "SEG: End Edit";
  editDialogOn = false;
  recoverDialog->hide();
  delete recoverDialog;
}

void EditSegmentorRecoverPlugin::Decorate(MeshModel &_m, GLArea * _gla) {
  qDebug() << "Decorate called";
  MeshlabDrawer* d = (MeshlabDrawer*)Segmentor::Instance()->getDrawer();
  if (d != NULL) {
	d->draw();
  }
}

void EditSegmentorRecoverPlugin::mousePressEvent(QMouseEvent *event, MeshModel &mm, GLArea *gla) {
  gla->suspendedEditor = true;
  QCoreApplication::sendEvent(gla, event);
  gla->suspendedEditor = false;
}

void EditSegmentorRecoverPlugin::mouseMoveEvent(QMouseEvent *event, MeshModel &mm, GLArea *gla) {
  gla->suspendedEditor = true;
  QCoreApplication::sendEvent(gla, event);
  gla->suspendedEditor = false;
}
void EditSegmentorRecoverPlugin::mouseReleaseEvent(QMouseEvent *event, MeshModel &mm, GLArea *gla) {
  gla->suspendedEditor = true;  
  QCoreApplication::sendEvent(gla, event);
  gla->suspendedEditor = false;
}
