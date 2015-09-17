#include <Qt>
#include <QSettings>
#include <QDir>
#include "edit_segmentor_recover.h"

using namespace std;
using namespace vcg;

EditSegmentorRecoverPlugin::EditSegmentorRecoverPlugin() {
	qFont.setFamily("Helvetica");
	qFont.setPixelSize(12);
	
	qDebug() << "Init Edit Segmentor";

	settingsFile = QDir::homePath()+ QDir::separator() + "segmentor.ini";
	loadSettings();
}

const QString EditSegmentorRecoverPlugin::Info() 
{
	return tr("Segmentor for Meshlab.");
}

bool EditSegmentorRecoverPlugin::StartEdit(MeshDocument &_md, GLArea *_gla ) {
	this->md = &_md;
	gla = _gla;
	QWidget* parent = gla->window();

	return false;
}
void EditSegmentorRecoverPlugin::loadSettings() {
	QSettings settings(settingsFile, QSettings::NativeFormat);
}

void EditSegmentorRecoverPlugin::saveSettings() {
	QSettings settings(settingsFile, QSettings::NativeFormat);
	settings.sync();
}
