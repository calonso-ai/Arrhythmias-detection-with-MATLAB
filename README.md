# Arrhythmias detection through Machine Learning and statistics

(See "Detección de arritmias mediante parámetros de calidad de señal y estadísticos" to complete project analysis)

The rate of false alarms on ICU is very high due to human factor error or bad qualiy of measurement instrument.

The objetive of the algorithm introduces on this document is, given a record set which generate alarms in some ICU, determinate  which one is produced by a truly arrhythmia and  which one are false alarms. In case of atrial fibrillation the only goal is to know which atrial fibrillation records contain arrhythmia.

This algorithm is characterized by phisiological signals features extraction and three clasifiers, one by arrhythmia; bradycardia, tachycardia, ventricular tahycardia to determine if exist false alarm (alarm produced by wrong measurement) or true alarm (alarm produced by truly alarm) and statisticals use to evaluate the existence of atrial fibrillation.

The range of detenction is 75-99%  for true alarms and  74-94% for false alarms. The best score is found for tachycardía and the worst for ventricular tachycardia.

(References in Wiki)
