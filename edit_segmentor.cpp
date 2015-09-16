#include <Qt>
#include <QSettings>
#include <QDir>
#include "edit_segmentor.h"

using namespace std;
using namespace vcg;

EditSegmentorPlugin::EditSegmentorPlugin() {
	qFont.setFamily("Helvetica");
	qFont.setPixelSize(12);
	
	qDebug() << "Init Edit Segmentor";

	settingsFile = QDir::homePath()+ QDir::separator() + "segmentor.ini";
	loadSettings();
}

const QString EditSegmentorPlugin::Info() 
{
	return tr("Segmentor for Meshlab.");
}

bool EditSegmentorPlugin::StartEdit(MeshDocument &_md, GLArea *_gla ) {
	this->md = &_md;
	gla = _gla;
	QWidget* parent = gla->window();
	return false;
}
void EditSegmentorPlugin::loadSettings() {
	QSettings settings(settingsFile, QSettings::NativeFormat);
}

void EditSegmentorPlugin::saveSettings() {
	QSettings settings(settingsFile, QSettings::NativeFormat);
	settings.sync();
}
