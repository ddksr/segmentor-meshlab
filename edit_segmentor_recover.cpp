#include <Qt>
#include <QSettings>

#include "edit_segmentor_recover.h"
#include "segRecoverDialog.h"

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
	delete recoverDialog;
  }

  recoverDialog = new segRecoverDialog(gla->window(), model, gla);
  recoverDialog->show();
  
  editDialogOn = true;
  
  return true;
}

void EditSegmentorRecoverPlugin::EndEdit(MeshModel &_md, GLArea *_gla) {
  qDebug() << "SEG: End Edit";
  editDialogOn = false;
  recoverDialog->hide();
  delete recoverDialog;
}
