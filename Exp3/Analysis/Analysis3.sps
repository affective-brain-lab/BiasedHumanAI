﻿* Encoding: UTF-8.
*AI images.
USE ALL.
COMPUTE filter_$=(control = 0).
VARIABLE LABELS filter_$ 'control = 0 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

*Generalized Linear Mixed Models. 
GENLINMIXED
  /DATA_STRUCTURE SUBJECTS=subject
  /FIELDS TARGET=choice TRIALS=NONE OFFSET=NONE
  /TARGET_OPTIONS DISTRIBUTION=MULTINOMIAL LINK=LOGIT
  /FIXED  EFFECTS=condition USE_INTERCEPT=TRUE
  /RANDOM EFFECTS=condition USE_INTERCEPT=TRUE SUBJECTS=subject COVARIANCE_TYPE=VARIANCE_COMPONENTS 
    SOLUTION=FALSE 
  /BUILD_OPTIONS TARGET_CATEGORY_ORDER=ASCENDING INPUTS_CATEGORY_ORDER=ASCENDING MAX_ITERATIONS=100 
    CONFIDENCE_LEVEL=95 DF_METHOD=SATTERTHWAITE COVB=MODEL PCONVERGE=0.000001(ABSOLUTE) SCORING=0 
    SINGULAR=0.000000000001
  /EMMEANS_OPTIONS SCALE=ORIGINAL PADJUST=LSD.

*Control.
USE ALL.
COMPUTE filter_$=(control = 1).
VARIABLE LABELS filter_$ 'control = 1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

*Generalized Linear Mixed Models. 
GENLINMIXED
  /DATA_STRUCTURE SUBJECTS=subject
  /FIELDS TARGET=choice TRIALS=NONE OFFSET=NONE
  /TARGET_OPTIONS DISTRIBUTION=MULTINOMIAL LINK=LOGIT
  /FIXED  EFFECTS=condition USE_INTERCEPT=TRUE
  /RANDOM EFFECTS=condition USE_INTERCEPT=TRUE SUBJECTS=subject COVARIANCE_TYPE=VARIANCE_COMPONENTS 
    SOLUTION=FALSE 
  /BUILD_OPTIONS TARGET_CATEGORY_ORDER=ASCENDING INPUTS_CATEGORY_ORDER=ASCENDING MAX_ITERATIONS=100 
    CONFIDENCE_LEVEL=95 DF_METHOD=SATTERTHWAITE COVB=MODEL PCONVERGE=0.000001(ABSOLUTE) SCORING=0 
    SINGULAR=0.000000000001
  /EMMEANS_OPTIONS SCALE=ORIGINAL PADJUST=LSD.