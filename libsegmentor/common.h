#ifndef LIBSEGMENTOR_COMMON_H
#define LIBSEGMENTOR_COMMON_H

struct RecoverySettings {
  float k2;
  float k3;
  float planarNormalCompatibility;
  int growingSteps;
  float regitonCompatibility;
  int maxNumOfSteps;
  float varianceOfNoise;
  float postProcessingMaxError;
  bool newDiscrepancy;
  bool usePostProcessing;
};

RecoverySettings* getDefaultRecoverySettings() {
  RecoverySettings conf = {
	1.1,
	1.1,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0
  };

  return &conf;
}

#endif
