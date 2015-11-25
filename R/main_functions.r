#'This function generates keyword and ignore lists based on the expansion of
#'concepts.
#'
#'@description Given any given concept_id or string of text this function
#'generates keyword and ignore lists based on the expansion of concepts (looking
#'at their synonyms).
#'
#'@param connection    The connection to the database server.
#'@param aphroditeConceptName  The string of text / concept name to use.
#'@param schema        The database schema being used.
#'@param dbms          The target DBMS for SQL to be rendered in.
#'
#'@details Takes the aphroditeConceptName looks for synonyms and builds a list
#'of related concepts using the vocabulary hierarchies
#'
#'@return A list with two elements: a list of positive keywords found
#'(keywordlist_ALL), and a list of ignore keywords (ignorelist_ALL)
#'
#' @examples \dontrun{
#'
#'wordLists <- buildKeywordList(conn, aphrodite_concept_name, cdmSchema, dbms)
#'
#' }
#'
#'@export
getBasicStatistics <- function(){
	cohort <- read.table("diabetesCohort")
	pid <- unique(cohort$patient)
	pidToKeep <- rep(0,length(pid))
	for(i in 1:length(pid)){
	  dat <- cohort[cohort$patient==pid[i],]
	  drugFreq <- as.data.frame(table(as.character(dat$drugid)))
	  if(length(drugFreq$Freq)>1) next
	  pidToKeep[i] <- 1
	}
	remove(drugFreq,dat)
	pid <- pid[which(pidToKeep==1)]
	remove(pidToKeep)
	#SubCohort of patient who remained on drug with what they started 
	cohortSameStartEnd <- data.frame()
	for(i in 1:length(pid)){
	  dat <- cohort[cohort$patient==pid[i],]
	  pidStartDate <- as.Date(dat$pidindexdate[1])
	  pidEndDate <- as.Date(dat$pidindexdate[length(dat$patient)])
	  pidStartDrug <- as.character(dat$drugid[1])
	  pidEndDrug <- as.character(dat$drugid[1])
	  pidGender <- dat$gender[1]
	  pidAge <- dat$age[1]
	  dat2 <- cbind(pid[i],as.character(pidStartDate),as.character(pidEndDate),pidStartDrug,pidEndDrug,pidGender,pidAge)
	  cohortSameStartEnd <- rbind(cohortSameStartEnd,dat2)
	}
	colnames(cohortSameStartEnd) <- c("patient","pidindexdate","pidenddate","pidstartdrug","pidenddrug","pidgender","pidage")

	#Case-2 Get the patient who started with one drug and move to second within at least 30 days ...
	cohort <- cohort[!(cohort$patient %in% pid),] #these are the patient who moved to other drugs
	remove(pid,dat2,pidStartDate,pidEndDate,pidStartDrug,pidEndDrug,pidGender,pidAge,dat)
	#SubCase-1 These patients should move to other drug not before 30 days !
	pid <- unique(cohort$patient)
	cohortStartSwitch <- data.frame()
	for(i in 1:length(pid)){
	  dat <- cohort[cohort$patient==pid[i],]
	  pidStartDate <- as.Date(dat$pidindexdate[1])
	  pidEndDate <- as.Date(dat$pidindexdate[length(dat$patient)])
	  pidStartDrug <- as.character(dat$drugid[1])
	  pidGender <- dat$gender[1]
	  pidAge <- dat$age[1]
	  for(j in 2:length(dat$patient)){
		x <- as.numeric(grep(paste("^",dat$drugid[j-1],"$",sep=""),dat$drugid[j]))
		if(length(x)!=1) break
	  }
	  if(j==2){
		pidMovedDate <- as.Date(dat$pidindexdate[j])
		pidMovedDrug <- as.character(dat$drugid[j])
	  }else
	  {
		pidMovedDate <- as.Date(dat$pidindexdate[j])
		pidMovedDrug <- as.character(dat$drugid[j])
	  }
	  remove(j)
	  if(as.numeric(pidMovedDate-pidStartDate)<30) next #Change here if you wanat to cut loose the constrain
	  dat2 <- cbind(pid[i],as.character(pidStartDate),as.character(pidMovedDate),as.character(pidEndDate),pidStartDrug,pidMovedDrug,pidGender,pidAge)
	  cohortStartSwitch <- rbind(cohortStartSwitch,dat2)
	}
	colnames(cohortStartSwitch) <- c("patient","pidindexdate","pidmoveddate","pidenddate","pidstartdrug","pidmoveddrug","pidgender","pidage")
	remove(pid,dat2,pidStartDate,pidEndDate,pidStartDrug,pidMovedDrug,pidMovedDate,x,i,cohort,pidGender,pidAge,dat)
	cohortSameStartEnd$drugPath <- paste(cohortSameStartEnd$pidstartdrug,cohortSameStartEnd$pidenddrug,sep="->")
	cohortStartSwitch$drugPath <- paste(cohortStartSwitch$pidstartdrug,cohortStartSwitch$pidmoveddrug,sep="->")
	pidDat1 <- subset(cohortSameStartEnd, select = c(patient,pidindexdate,pidgender,pidage,drugPath))
	pidDat2 <- subset(cohortStartSwitch, select = c(patient,pidindexdate,pidgender,pidage,drugPath))
	pidDat <- rbind(pidDat1,pidDat2)
	remove(pidDat1,pidDat2)

	# Get total number of patients and their age distribution in the cohort overall and then for each treatment pathway (kind of matrix plot in ggplot)
	# 1) Getting age for all
	pidDat$pidage <- as.numeric(as.character(pidDat$pidage))
	allPidAge <- ggplot(pidDat, aes(x=pidage)) + geom_histogram(color="white",fill="darkblue",binwidth=4, alpha=0.85) + theme_bw() + theme(text = element_text(size=16), axis.text.x = element_text(angle=0, vjust=1)) + labs(x = "Age", y = "Patients",title = "Age Distribution of all Patients") 

	# 2) Age dist of patient who started with one and remained on one
	cohortSameStartEnd$pidage <- as.numeric(as.character(cohortSameStartEnd$pidage))
	allPidAgeSameStartEnd <- ggplot(cohortSameStartEnd, aes(x=pidage)) + geom_histogram(color="white",fill="darkblue",binwidth=4, alpha=0.85) + theme_bw() + theme(text = element_text(size=16), axis.text.x = element_text(angle=0, vjust=1)) + labs(x = "Age", y = "Number of Patients", title = "Age Dist - Patients with Same Start End Drug")
	# 3) Age dist patient who started with one and moved to other
	cohortStartSwitch$pidage <- as.numeric(as.character(cohortStartSwitch$pidage))
	allPidAgeStartSwitch <- ggplot(cohortStartSwitch, aes(x=pidage)) + geom_histogram(color="white",fill="darkblue",binwidth=4, alpha=0.85) + theme_bw() + theme(text = element_text(size=16), axis.text.x = element_text(angle=0, vjust=1)) + labs(x = "Age", y = "Number of Patients", title = "Age Dist - Patients with Start and Moved to Other Drug")
	# Age distribution according to treatment pathways
	drugPair <- as.data.frame(table(pidDat$drugPath))
	colnames(drugPair) <- c("drugPath","Freq")
	drugPair <- drugPair[order(-drugPair$Freq),]
	drugPair$drugPath <- as.character(drugPair$drugPath)
	#Should have atleast 50 patients for histogram
	drugPair <- filter(drugPair,Freq>=50)
	drugAge <- data.frame()
	for(i in 1:nrow(drugPair)){
	  dat <- pidDat[pidDat$drugPath==drugPair$drugPath[i],]
	  drugAge <- rbind(drugAge,dat)
	}
	drugAge$pidage <- as.numeric(as.character(drugAge$pidage))
	drugAge <- subset(drugAge,select=c(drugPath,pidage))
	drugAgePlot <- ggplot(drugAge) + geom_density(aes(x = pidage,fill = (drugPath)),alpha=0.4) + labs(x = NULL) + theme_bw() + theme(text = element_text(size=16), axis.text.x = element_text(angle=0, vjust=1)) + labs(x = "Age", y = "Number of Patients",title = "Age Density - All Pair Treatment Pathways")

	#Aggregate the data to plot it into single pdf file
	dataToPlot <- list(allPidAge,allPidAgeSameStartEnd,allPidAgeStartSwitch,drugAgePlot)
	pdf("ageStatPatients.pdf")
	for(i in 1:length(dataToPlot)){
	  print(dataToPlot[[i]])
	}
	dev.off()
}
#'This function generates keyword and ignore lists based on the expansion of
#'concepts.
#'
#'@description Given any given concept_id or string of text this function
#'generates keyword and ignore lists based on the expansion of concepts (looking
#'at their synonyms).
#'
#'@param connection    The connection to the database server.
#'@param aphroditeConceptName  The string of text / concept name to use.
#'@param schema        The database schema being used.
#'@param dbms          The target DBMS for SQL to be rendered in.
#'
#'@details Takes the aphroditeConceptName looks for synonyms and builds a list
#'of related concepts using the vocabulary hierarchies
#'
#'@return A list with two elements: a list of positive keywords found
#'(keywordlist_ALL), and a list of ignore keywords (ignorelist_ALL)
#'
#' @examples \dontrun{
#'
#'wordLists <- buildKeywordList(conn, aphrodite_concept_name, cdmSchema, dbms)
#'
#' }
#'
#'@export
getCohort <- function(conn){
  conn <- connect(connectionDetails)
  t2dDrug <- c("173","10633","2404","4821","217360","4815","25789","73044","274332","6809","84108","33738","72610","16681","30009","593411","60548")
  t2dDrug <- paste(t2dDrug,sep=",",collapse=",")
  t2dDrug <- gsub(",","','",t2dDrug)
  queryMe <- paste("SELECT concept_ID as drugid, concept_name as drugname FROM @cdmSchema.concept WHERE concept_code IN ('",t2dDrug,"') AND vocabulary_id = 'RxNorm' AND domain_ID = 'Drug';", sep="")
  drugInfo <- executeSQL_ro(conn, cdmSchema, queryMe, dbms)
  drugInfo$drugname <- toupper(drugInfo$drugname)
  drugInfo$shortname <- substring(drugInfo$drugname,1,4)
  remove(queryMe)
#  dbDisconnect(conn)
  drugConcept <- paste(drugInfo$drugid,sep=",",collapse=",")
  drugConcept <- gsub(",","','",drugConcept)
  queryMe <- paste("SELECT person_id as patient, drug_concept_id as drugid, drug_exposure_start_date as pidindexdate FROM @cdmSchema.drug_exposure WHERE drug_concept_id in ('",drugConcept,"')",sep="")
  pidDrug <- executeSQL_ro(conn, cdmSchema, queryMe, dbms)
#  dbDisconnect(conn)
  #Geeting first prescription of each of these patients
  pid <- unique(pidDrug$patient)
  patientDrug <- data.frame()
  for(i in 1:length(pid)){
    dat <- pidDrug[pidDrug$patient==pid[i],]
    dat <- dat[order(dat$pidindexdate),]
    dat <- aggregate(drugid~patient+pidindexdate,toString,data=dat)
    dat$drugid <- gsub(" ","",dat$drugid)
    dat2 <- sapply(dat$drugid, function(x) {paste(unique(sort(as.numeric(unlist(strsplit(x,","))))),sep="-",collapse="-")})
    dat2 <- as.data.frame(dat2)
    rownames(dat2) <- NULL
    dat$drugid <- dat2$dat2
    patientDrug <- rbind(patientDrug,dat)
  }
  #Replacing drugid with drugname
  for(i in 1:length(drugInfo$drugid)){
    patientDrug$drugid <- gsub(drugInfo$drugid[i],drugInfo$shortname[i],patientDrug$drugid)
  }
  remove(pidDrug,dat,dat2)
  #Removing patients who do not have follow-up data
  pidToKeep <- as.data.frame(table(patientDrug$patient))
  x <- which(pidToKeep$Freq>1)
  pidToKeep <- pidToKeep[x,]
  remove(x)
  x <- which(patientDrug$patient %in% pidToKeep$Var1)
  patientDrug <- patientDrug[x,]
  remove(x,pidToKeep)
  #Remove patients who do not have at least 90 days of data prior to their index date
  #conn <- connect(connectionDetails)
  pid <- paste(unique(patientDrug$patient),sep=",",collapse=",")
  pid <- gsub(",","','",pid)
  queryMe <- paste("SELECT person_id as patient, observation_period_start_date as pidstart, observation_period_end_date as pidend FROM @cdmSchema.observation_period WHERE person_id IN ('",pid,"')",sep="")
  pidObservation <- executeSQL_ro(conn, cdmSchema, queryMe, dbms)
  #dbDisconnect(conn)
  remove(pid,queryMe)
  pidDat <- merge(patientDrug,pidObservation,by="patient")
  pid <- unique(pidDat$patient)
  pidDat2 <- data.frame()
  for (i in 1:length(pid)){
    dat <- pidDat[pidDat$patient==pid[i],]
    dat <- dat[order(dat$pidindexdate),]
    dat <- dat[1,]
    dat$days <- as.numeric(as.Date(dat$pidindexdate) - as.Date(dat$pidstart))
    pidDat2 <- rbind(pidDat2,dat)
  }
  remove(dat,pidDat,i)
  x <- which(pidDat2$days>=90)
  pidDat2 <- pidDat2[x,]
  remove(x)
  x <- which(patientDrug$patient %in% pidDat2$patient)
  patientDrug <- patientDrug[x,]
  remove(x,pidDat2)
  patientDrug <- merge(patientDrug,pidObservation,by="patient")
  # Removing patient who were treated with any of the Type-1 diabetes drug before their index time
  #conn <- connect(connectionDetails)
  t1dDrug <- c("139825","274783","314684","352385","400008","51428","5856","86009","139953")
  t1dDrug <- paste(t1dDrug,sep=",",collapse=",")
  t1dDrug <- gsub(",","','",t1dDrug)
  queryMe <- paste("SELECT concept_id as drugid, concept_name as drugname from @cdmSchema.concept WHERE concept_code IN ('",t1dDrug,"') and Vocabulary_id = 'RxNorm' and domain_id = 'Drug';", sep="")
  drugInfoT1D <- executeSQL_ro(conn, cdmSchema, queryMe, dbms)
  remove(queryMe)
  drugConceptT1D <- paste(drugInfoT1D$drugid,sep=",",collapse=",")
  drugConceptT1D <- gsub(",","','",drugConceptT1D)
  queryMe <- paste("select person_id as patient, drug_concept_id as drugid, drug_exposure_start_date as pidindexdate FROM @cdmSchema.drug_exposure WHERE drug_concept_id IN ('",drugConceptT1D,"');",sep="")
  pidDrugT1D <- executeSQL_ro(conn, cdmSchema, queryMe, dbms)
  #dbDisconnect(conn)
  pidDrugT1D <- merge(pidDrugT1D,drugInfoT1D,by="drugid")
  pid <- unique(patientDrug$patient)
  pidDat <- data.frame()
  for(i in 1:length(pid)){
    dat <- patientDrug[patientDrug$patient==pid[i],]
    dat <- dat[1,]
    pid2 <- dat$patient
    dat2 <- pidDrugT1D[pidDrugT1D$patient==pid2,]
    if(dim(dat2)[1]==0) next
    dat2 <- subset(dat2,as.Date(dat2$pidindexdate) < as.Date(dat$pidindexdate))
    pidDat <- rbind(pidDat,dat2)
  }
  pid <- unique(pidDat$patient)
  x <- which(patientDrug$patient %in% pid)
  patientDrug <- patientDrug[-x,]
  remove(pidDat,pidDrugT1D,pidObservation)
  #Getting demographic information for these patients
  #conn <- connect(connectionDetails)
  pid <- as.character(unique(patientDrug$patient))
  pid <- paste(pid,sep=",",collapse=",")
  pid <- gsub(",","','",pid)
  queryMe <- paste("SELECT person_id as patient, gender_concept_id as gender, year_of_birth as birthyear, race_concept_id as race FROM @cdmSchema.person WHERE person_id IN ('",pid,"');",sep="")
  pidDem <- executeSQL_ro(conn, cdmSchema, queryMe, dbms)
  #dbDisconnect(conn)
  patientDrug <- merge(patientDrug,pidDem,by="patient")
  dat <- patientDrug
  dat$age <- as.numeric(as.numeric((format(as.Date(patientDrug$pidstart),format="%Y")))-as.numeric(as.character(patientDrug$birthyear)))
  patientDrug <- dat
  diabetesCohort <- patientDrug
  write.table(diabetesCohort,file="diabetesCohort")
  #Getting the measurement data (HbA1c)
  #conn <- connect(connectionDetails)
  queryMe <- paste("SELECT person_id as patient, measurement_date as measureddate, value_as_number as labval FROM @cdmSchema.measurement WHERE measurement_concept_id = 40789263;",sep="")
  measurements <- executeSQL_ro(conn, cdmSchema, queryMe, dbms)
  write.table(measurements,file="measurements")
  #dbDisconnect(conn)
}
#'This function generates keyword and ignore lists based on the expansion of
#'concepts.
#'
#'@description Given any given concept_id or string of text this function
#'generates keyword and ignore lists based on the expansion of concepts (looking
#'at their synonyms).
#'
#'@param connection    The connection to the database server.
#'@param aphroditeConceptName  The string of text / concept name to use.
#'@param schema        The database schema being used.
#'@param dbms          The target DBMS for SQL to be rendered in.
#'
#'@details Takes the aphroditeConceptName looks for synonyms and builds a list
#'of related concepts using the vocabulary hierarchies
#'
#'@return A list with two elements: a list of positive keywords found
#'(keywordlist_ALL), and a list of ignore keywords (ignorelist_ALL)
#'
#' @examples \dontrun{
#'
#'wordLists <- buildKeywordList(conn, aphrodite_concept_name, cdmSchema, dbms)
#'
#' }
#'
#'@export
getCombinedFeatures <- function(conn){
	message("Building Patient Feature Matrix")
	conn <- connect(connectionDetails)
	cohort <- read.table("diabetesCohort")
	pid <- unique(cohort$patient)
	pidToKeep <- rep(0,length(pid))
	for(i in 1:length(pid)){
	  dat <- cohort[cohort$patient==pid[i],]
	  drugFreq <- as.data.frame(table(as.character(dat$drugid)))
	  if(length(drugFreq$Freq)>1) next
	  pidToKeep[i] <- 1
	}
	remove(drugFreq,dat)
	pid <- pid[which(pidToKeep==1)]
	remove(pidToKeep)
	#SubCohort of patient who remained on drug with what they started with 
	cohortSameStartEnd <- data.frame()
	for(i in 1:length(pid)){
	  dat <- cohort[cohort$patient==pid[i],]
	  pidStartDate <- as.Date(dat$pidindexdate[1])
	  pidEndDate <- as.Date(dat$pidindexdate[length(dat$patient)])
	  pidStartDrug <- as.character(dat$drugid[1])
	  pidEndDrug <- as.character(dat$drugid[1])
	  pidGender <- dat$gender[1]
	  pidAge <- dat$age[1]
	  dat2 <- cbind(pid[i],as.character(pidStartDate),as.character(pidEndDate),pidStartDrug,pidEndDrug,pidGender,pidAge)
	  cohortSameStartEnd <- rbind(cohortSameStartEnd,dat2)
	}
	colnames(cohortSameStartEnd) <- c("patient","pidindexdate","pidenddate","pidstartdrug","pidenddrug","pidgender","pidage")
	#Case-2 Get the patient who started with one drug and moved to second within at least 30 days ...
	cohort <- cohort[!(cohort$patient %in% pid),]
	remove(pid,dat2,pidStartDate,pidEndDate,pidStartDrug,pidEndDrug,pidGender,pidAge,dat)
	pid <- unique(cohort$patient)
	cohortStartSwitch <- data.frame()
	for(i in 1:length(pid)){
	  dat <- cohort[cohort$patient==pid[i],]
	  pidStartDate <- as.Date(dat$pidindexdate[1])
	  pidEndDate <- as.Date(dat$pidindexdate[length(dat$patient)])
	  pidStartDrug <- as.character(dat$drugid[1])
	  pidGender <- dat$gender[1]
	  pidAge <- dat$age[1]
	  for(j in 2:length(dat$patient)){
		x <- as.numeric(grep(paste("^",dat$drugid[j-1],"$",sep=""),dat$drugid[j]))
		if(length(x)!=1) break
	  }
	  if(j==2){
		pidMovedDate <- as.Date(dat$pidindexdate[j])
		pidMovedDrug <- as.character(dat$drugid[j])
	  }else
	  {
		pidMovedDate <- as.Date(dat$pidindexdate[j])
		pidMovedDrug <- as.character(dat$drugid[j])
	  }
	  remove(j)
	  if(as.numeric(pidMovedDate-pidStartDate)<30) next
	  dat2 <- cbind(pid[i],as.character(pidStartDate),as.character(pidMovedDate),as.character(pidEndDate),pidStartDrug,pidMovedDrug,pidGender,pidAge)
	  cohortStartSwitch <- rbind(cohortStartSwitch,dat2)
	}
	colnames(cohortStartSwitch) <- c("patient","pidindexdate","pidmoveddate","pidenddate","pidstartdrug","pidmoveddrug","pidgender","pidage")
	remove(pid,dat2,pidStartDate,pidEndDate,pidStartDrug,pidMovedDrug,pidMovedDate,x,i,cohort,pidGender,pidAge,dat)
	cohortSameStartEnd$drugPath <- paste(cohortSameStartEnd$pidstartdrug,cohortSameStartEnd$pidenddrug,sep="->")
	cohortStartSwitch$drugPath <- paste(cohortStartSwitch$pidstartdrug,cohortStartSwitch$pidmoveddrug,sep="->")
	pidDat1 <- subset(cohortSameStartEnd, select = c(patient,pidindexdate,pidgender,pidage,drugPath))
	pidDat2 <- subset(cohortStartSwitch, select = c(patient,pidindexdate,pidgender,pidage,drugPath))
	pidDat <- rbind(pidDat1,pidDat2)
	remove(cohortSameStartEnd,cohortStartSwitch,pidDat1,pidDat2)

	# Geeting Drug Feature
	drugPair <- as.data.frame(table(pidDat$drugPath))
	drugPair <- drugPair[order(-drugPair$Freq),]
	drugPair <- filter(drugPair,Freq>=50) #Keeping pairs with more than or equal to 50 patients
	colnames(drugPair) <- c("drugPath")

	# Generate the list of matrices - pair wise that we can use for random forest
	pair <- t(combn(length(drugPair$drugPath),2))
	pidDrugBeforeIndex <- list()
	for(i in 1:nrow(pair)){
	  dat1 <- pidDat[pidDat$drugPath==drugPair$drugPath[pair[i,1]],]
	  dat2 <- pidDat[pidDat$drugPath==drugPair$drugPath[pair[i,2]],]
	  pidDatFinal <- rbind(dat1,dat2)
	  pidDrugBefIndex <- data.frame()
	  pid <- as.numeric(as.character(pidDatFinal$patient))
	  for (j in 1:length(pid)){
		dat <- pidDatFinal[pidDatFinal$patient==pid[j],]
		indexDate <- as.Date(dat$pidindexdate)
		pid2 <- dat$patient
		queryMe <- paste("SELECT person_id AS patient, drug_concept_id AS drugid, drug_exposure_start_date AS datebeforeindex  FROM @cdmSchema.drug_exposure WHERE person_id = ",pid2," AND drug_exposure_start_date < '",indexDate,"' ORDER BY datebeforeindex;", sep="")
		dat2 <- executeSQL_ro(conn, cdmSchema, queryMe, dbms)
		pidDrugBefIndex <- rbind(pidDrugBefIndex,dat2)
	  }
	  pidDrugBefIndex <- pidDrugBefIndex[,-3]
	  remove(dat1,dat2,pidDatFinal,pid,dat,indexDate,pid2,queryMe,dat2)
	  pidDrugBeforeIndex[[i]] <- pidDrugBefIndex
	}
	#dbDisconnect(conn)
	pidDrug <- list()
	for(i in 1:length(pidDrugBeforeIndex)){
	  dat <- pidDrugBeforeIndex[[i]]
	  dat$drugid <- paste("Rx",dat$drugid,sep="")
	  dat <- dat[!duplicated(dat),]
	  pidDrug[[i]] <- dat
	}
	# Getting Dx Features
	pidConditionBeforeIndex <- list()
	conn <- connect(connectionDetails)
	for(i in 1:nrow(pair)){
	  dat1 <- pidDat[pidDat$drugPath==drugPair$drugPath[pair[i,1]],]
	  dat2 <- pidDat[pidDat$drugPath==drugPair$drugPath[pair[i,2]],]
	  pidDatFinal <- rbind(dat1,dat2)
	  pidCondBefIndex <- data.frame()
	  pid <- as.numeric(as.character(pidDatFinal$patient))
	  for (j in 1:length(pid)){
		dat <- pidDatFinal[pidDatFinal$patient==pid[j],]
		indexDate <- as.Date(dat$pidindexdate)
		pid2 <- dat$patient
		queryMe <- paste("SELECT person_id AS patient, condition_concept_id AS conditionid, condition_start_date AS datebeforeindex  FROM @cdmSchema.condition_occurrence WHERE person_id = ",pid2," AND condition_start_date < '",indexDate,"' ORDER BY datebeforeindex;", sep="")
		dat2 <- executeSQL_ro(conn, cdmSchema, queryMe, dbms)
		pidCondBefIndex <- rbind(pidCondBefIndex,dat2)
	  }
	  pidCondBefIndex <- pidCondBefIndex[,-3]
	  remove(dat1,dat2,pidDatFinal,pid,dat,indexDate,pid2,queryMe,dat2)
	  pidConditionBeforeIndex[[i]] <- pidCondBefIndex
	}
	#dbDisconnect(conn)  
	pidCondition <- list()
	for(i in 1:length(pidConditionBeforeIndex)){
	  dat <- pidConditionBeforeIndex[[i]]
	  dat$conditionid <- paste("Dx",dat$conditionid,sep="")
	  dat <- dat[!duplicated(dat),]
	  pidCondition[[i]] <- dat
	}
	pidFeatureMatDat <- list()
	for(i in 1:length(pidCondition)){
	  dat1 <- pidDrug[[i]]
	  dat2 <- pidCondition[[i]]
	  colnames(dat1) <- c("patient","feature")
	  colnames(dat2) <- c("patient","feature")
	  dat <- rbind(dat1,dat2)
	  dat$val <- rep(1,length(dat$patient))
	  datMat <- acast(dat,patient~feature,value.var="val")
	  datMat[is.na(datMat)] <- 0
	  datMat <- cbind(rownames(datMat),datMat)
	  colnames(datMat)[1] <- c("patient")
	  datFeatureMat <- merge(datMat,pidDat,by="patient")
	  pidFeatureMatDat[[i]] <- datFeatureMat
	  remove(dat,datMat,datFeatureMat)
	}
	#save(pidFeatureMatDat,file="pidFeatureMatDat") # Can save it here to avoide building the same matrix from SQL every time you run the code.
	remove(dat1,dat2,drugPair,pair,pidCondBefIndex,pidDat,pidDrugBefIndex,i,j,pidCondition,pidConditionBeforeIndex,pidDrug,pidDrugBeforeIndex)
	# Building classifier
	registerDoMC(4) # Can change here based on the CPU config.
	curClassifier <- list()
	predRocFinal <- list()
	myroc <- list()
	curConfusionMatrix <- list()
	imp <- list()
	drugPath <- list()
	threshold <- list()
	posClass <- list()
	negClass <- list()
	for(i in 1:length(pidFeatureMatDat)){
	  dat <- pidFeatureMatDat[[i]]
	  dat$drugPath <- gsub("->","",dat$drugPath)
	  drugs <- unique(dat$drugPath)
	  d1 <- as.character(drugs[1])
	  d2 <- as.character(drugs[2])
	  dat2 <- dat[dat$drugPath==d1,] #1
	  dat3 <- dat[dat$drugPath==d2,] #2
	  x <- c(length(dat2$patient),length(dat3$patient))
	  bigIndex <- which(x==max(x))
	  smallIndex <- which(x==min(x))
	  set.seed(440)
	  if(bigIndex==1){
		dat2 <- sample_n(dat2,length(dat3$drugPath),replace=TRUE)
	  }else
	  {
		if(bigIndex==2){
		  dat3 <- sample_n(dat3,length(dat2$drugPath),replace=TRUE)
		}
	  }
	  dat4 <- rbind(dat2,dat3)
	  dat5 <- subset(dat4,select = c(pidage,pidgender,drugPath))
	  dat4 <- subset(dat4, select = -c(pidage,pidgender,drugPath,patient,pidindexdate))  
	  col <- c(1:ncol(dat4))  
	  dat4[,col]<-lapply(col, function(x) as.numeric(as.character(dat4[,x])))
	  dat4[dat4>=1]<-1
	  colmnSum <- as.numeric(colSums(dat4))
	  xx <- which(colmnSum>=10)
	  dat4 <- dat4[,c(xx)]
	  datFinal <- cbind(dat4,dat5)
	  xxx <- which(colnames(datFinal)=="Dx201826")
	  datFinal <- datFinal[,-xxx]
	  datFinal$pidage <- as.numeric(as.character(datFinal$pidage))
	  datFinal$pidgender <- factor(datFinal$pidgender)
	  datFinal$drugPath <- factor(datFinal$drugPath)
	  datFinal <- as.data.frame(datFinal)
	  datPat <- createDataPartition(datFinal$drugPath, p = .80, list = FALSE)
	  training <- datFinal[datPat,]
	  testing  <- datFinal[-datPat,]
	  predLabel <- subset(training, select = -c(drugPath))
	  predLabel <- colnames(predLabel)
	  cvCtrl = trainControl(method = "repeatedcv",number = 10, repeats = 5, classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = TRUE)
	  newGrid = expand.grid(mtry = c(2,4,8,15))
	  posClass[[i]] <- d1
	  negClass[[i]] <- d2
	  classifierRF = train(drugPath ~ ., data = training, trControl = cvCtrl, method = "rf", metric="ROC", tuneGrid = newGrid)
	  curClassifier[[i]] <- classifierRF
	  predRoc = predict(classifierRF, testing, type = "prob")
	  predRocFinal[[i]] <- predRoc
	  myroc[[i]] = pROC::roc(testing$drugPath, as.vector(predRoc[,2]))
	  threshold[[i]] = coords(myroc[[i]],x="best",best.method = "closest.topleft")[[1]]
	  predCut = factor( ifelse(predRoc[, d1] > threshold[[i]], d1, d2) )
	  curConfusionMatrix[[i]] = confusionMatrix(predCut, testing$drugPath, positive = d1)
	  imp[[i]] <- varImp(classifierRF, scale = FALSE) 
	  drugPath[[i]] <- drugs
	  remove(dat,dat2,d1,d2,dat3,dat4,dat5,datFinal,training,testing,x,xx,xxx,classifierRF,predRoc,predCut,datPat,predLabel,cvCtrl,newGrid)
	}

	# Performing Lasso
	modelLasso <- list()
	myrocLasso <- list()
	curConfusionMatrixLasso <- list()
	impLasso <- list()
	predRocLassoFinal <- list()
	thresholdLasso <- list()
	for(i in 1:length(pidFeatureMatDat)){
	  dat <- pidFeatureMatDat[[i]]
	  dat$drugPath <- gsub("->","",dat$drugPath)
	  drugs <- unique(dat$drugPath)
	  d1 <- as.character(drugs[1])
	  d2 <- as.character(drugs[2])
	  dat2 <- dat[dat$drugPath==d1,] #1
	  dat3 <- dat[dat$drugPath==d2,] #2
	  x <- c(length(dat2$patient),length(dat3$patient))
	  bigIndex <- which(x==max(x))
	  smallIndex <- which(x==min(x))
	  set.seed(440)
	  if(bigIndex==1){
		dat2 <- sample_n(dat2,length(dat3$drugPath),replace=TRUE)
	  }else
	  {
		if(bigIndex==2){
		  dat3 <- sample_n(dat3,length(dat2$drugPath),replace=TRUE)
		}
	  }
	  dat4 <- rbind(dat2,dat3)
	  dat5 <- subset(dat4,select = c(pidage,pidgender,drugPath))
	  dat4 <- subset(dat4, select = -c(pidage,pidgender,drugPath,patient,pidindexdate))  
	  col <- c(1:ncol(dat4))  
	  dat4[,col]<-lapply(col, function(x) as.numeric(as.character(dat4[,x])))
	  dat4[dat4>=1]<-1
	  colmnSum <- as.numeric(colSums(dat4))
	  xx <- which(colmnSum>=10)
	  dat4 <- dat4[,c(xx)]
	  datFinal <- cbind(dat4,dat5)
	  xxx <- which(colnames(datFinal)=="Dx201826")
	  datFinal <- datFinal[,-xxx]
	  datFinal$pidage <- as.numeric(as.character(datFinal$pidage))
	  datFinal$pidgender <- factor(datFinal$pidgender)
	  datFinal$drugPath <- factor(datFinal$drugPath)
	  datFinal <- as.data.frame(datFinal)
	  datPat <- createDataPartition(datFinal$drugPath, p = .80, list = FALSE)
	  training <- datFinal[datPat,]
	  testing  <- datFinal[-datPat,]
	  predLabel <- subset(training, select = -c(drugPath))
	  predLabel <- colnames(predLabel)
	  cvCtrl = trainControl(method = "repeatedcv",number = 10, repeats = 5, classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = TRUE)
	  gridLasso <- expand.grid(alpha = (1:10) * 0.1, lambda = (1:10) * 0.1)
	  modeLasso <- train(drugPath ~ ., data = training, method = "glmnet", family = "binomial", trControl = cvCtrl, metric = "ROC", tuneGrid = gridLasso)
	  modelLasso[[i]] <- modeLasso
	  predRocLasso <- predict(modeLasso, testing, type = "prob")
	  predRocLassoFinal[[i]] <- predRocLasso
	  myrocLasso[[i]] = pROC::roc(testing$drugPath, as.vector(predRocLasso[,2]))
	  thresholdLasso[[i]] = coords(myrocLasso[[i]],x="best",best.method = "closest.topleft")[[1]]
	  predCutLasso = factor( ifelse(predRocLasso[, d1] > thresholdLasso[[i]], d1, d2) )
	  curConfusionMatrixLasso[[i]] = confusionMatrix(predCutLasso, testing$drugPath, positive = d1)
	  impLasso[[i]] <- varImp(modeLasso, scale = FALSE) 
	  remove(dat,dat2,d1,d2,dat3,dat4,dat5,datFinal,training,testing,x,xx,xxx,modeLasso,predCutLasso,datPat,predLabel,cvCtrl,gridLasso)
	}
	modStatRf <- data.frame()
	for(i in 1:length(curConfusionMatrix)){
	  datClassifier <- curClassifier[[i]]
	  datMat <- curConfusionMatrix[[i]]
	  datROC <- myroc[[i]]
	  dat <- cbind(posClass[[i]],negClass[[i]],datROC$auc[[1]],datMat$byClass[[1]],datMat$byClass[[2]],datMat$byClass[[3]])
	  modStatRf <- rbind(modStatRf,dat)
	}
	colnames(modStatRf) <- c("classOne","classTwo","AUC","sensitivity","spectificity","ppv")
	modStatRf$AUC <- format(round(as.numeric(as.character(modStatRf$AUC)),2))
	modStatRf$sensitivity <- format(round(as.numeric(as.character(modStatRf$sensitivity)),2))
	modStatRf$spectificity <- format(round(as.numeric(as.character(modStatRf$spectificity)),2))
	modStatRf$ppv <- format(round(as.numeric(as.character(modStatRf$ppv)),2))
	# For Lasso - Use this in case want to report for Lasso
	#modStatLasso <- data.frame()
	#for(i in 1:length(curConfusionMatrixLasso)){
	#  datClassifier <- modelLasso[[i]]
	#  datMat <- curConfusionMatrixLasso[[i]]
	#  datROC <- myrocLasso[[i]]
	#  dat <- cbind(posClass[[i]],negClass[[i]],datROC$auc[[1]],datMat$byClass[[1]],datMat$byClass[[2]],datMat$byClass[[3]])
	#  modStatLasso <- rbind(modStatLasso,dat)
	#}
	#colnames(modStatLasso) <- c("classOne","classTwo","AUC","sensitivity","spectificity","ppv")
	#modStatLasso$AUC <- format(round(as.numeric(as.character(modStatLasso$AUC)),2))
	#modStatLasso$sensitivity <- format(round(as.numeric(as.character(modStatLasso$sensitivity)),2))
	#modStatLasso$spectificity <- format(round(as.numeric(as.character(modStatLasso$spectificity)),2))
	#modStatLasso$ppv <- format(round(as.numeric(as.character(modStatLasso$ppv)),2))

	#Ploting ROC curves
	pdf("modelRocCurves.pdf")
	for(i in 1:length(myroc)){
	  plot.roc(smooth(myroc[[i]]),col="black")
	  plot.roc(smooth(myrocLasso[[i]]),col="red",add=TRUE)
	  legend("bottomright", legend=c("Random Forest", "Lasso"),col=c(par("fg"), "red"), lwd=2)
	  title(main = paste(myroc[[i]]$levels[2],myroc[[i]]$levels[1],sep=" & "),line=2.5)
	}
	dev.off()
	pdf("modelRfStatistics.pdf")
	grid.table(modStatRf)
	dev.off()
	# Geeting important feature matrix plot. This is only for the Random forest model. Do same if you want for Lasso
	# Only reporting the top 10 features for each drug pair.
	impFeatureClass <- list()
	conn <- connect(connectionDetails)
	for(i in 1:length(posClass)){
	  dat <- imp[[i]]
	  impFeature <- data.frame(dat$importance)
	  impFeature <- cbind(rownames(impFeature),impFeature)
	  colnames(impFeature) <- c("feature","importance")
	  impFeature <- impFeature[order(-impFeature$importance),]
	  impFeature$feature <- gsub("pidage","4265453",impFeature$feature)
	  impFeature$feature <- gsub("pidgender8532","8532",impFeature$feature)
	  impFeature$feature <- gsub("pidgender8507","8507",impFeature$feature)
	  impFeature <- impFeature[1:15,] # selecting top 15 only
	  impFeature$feature <- gsub("Dx","",impFeature$feature)
	  impFeature$feature <- gsub("Rx","",impFeature$feature)
	  impFeature <- cbind(rownames(impFeature),impFeature)
	  colnames(impFeature) <- c("featureDxRx","feature","importance")
	  allFeatures <- paste(impFeature$feature,sep=",",collapse=",")
	  allFeatures <- gsub(",","','",allFeatures)
	  queryMe <- paste("SELECT concept_id AS feature, concept_name AS featurename FROM @cdmSchema.concept WHERE concept_id IN ('",allFeatures,"');", sep="")
	  dat2 <- executeSQL_ro(conn, cdmSchema, queryMe, dbms)
	  impFeature <- merge(impFeature,dat2,by="feature")
	  impFeature <- impFeature[order(-impFeature$importance),]
	  impFeature$posclass <- posClass[[i]]
	  impFeature$negclass <- negClass[[i]]
	  impFeatureClass[[i]] <- impFeature
	  remove(dat2,queryMe)
	}
	dbDisconnect(conn) 
	remove(impFeature,allFeatures)
	# Generating matrix
	impFeatureMat <- data.frame()
	for(i in 1:length(impFeatureClass)){
	  dat <- impFeatureClass[[i]]
	  dat <- subset(dat,select=c(posclass,negclass,featurename))
	  impFeatureMat <- rbind(impFeatureMat,dat)
	}
	remove(dat)
	impFeatureMat$drugPath <- paste(impFeatureMat$posclass,impFeatureMat$negclass,sep="->")
	impFeatureMat <- subset(impFeatureMat,select=c(drugPath,featurename))
	uniqueFeature <- as.character(unique(impFeatureMat$featurename))
	uniqueDrugPath <- as.character(unique(impFeatureMat$drugPath))
	matToPlot <- data.frame()
	for(i in 1:length(uniqueDrugPath)){
	  featureVec <- rep(0,length(uniqueFeature))
	  dat <- impFeatureMat[impFeatureMat$drugPath==uniqueDrugPath[[i]],]
	  dat$featurename <- as.character(dat$featurename)
	  x <- which(uniqueFeature %in% dat$featurename)
	  featureVec[x] <- 1
	  dat2 <- cbind(uniqueDrugPath[[i]],uniqueFeature,featureVec)
	  dat2 <- as.data.frame(dat2)
	 # colnames(dat2) <- as.character(c(1:(1+length(uniqueFeature))))
	  matToPlot <- rbind(matToPlot,dat2)
	 remove(featureVec)
	}
	colnames(matToPlot) <- c("drugPath","feature","value")
	p <- ggplot(matToPlot, aes(feature, drugPath,fill=factor(value))) + geom_tile() + scale_fill_manual(values=c("0"="white", "1"="maroon")) + theme_bw() + theme(text = element_text(size=6), axis.text.x = element_text(angle=90, hjust=1))
	pdf("finalFeatureMatrix.pdf")
	print(p)
	dev.off()
}
#'This function generates keyword and ignore lists based on the expansion of
#'concepts.
#'
#'@description Given any given concept_id or string of text this function
#'generates keyword and ignore lists based on the expansion of concepts (looking
#'at their synonyms).
#'
#'@param connection    The connection to the database server.
#'@param aphroditeConceptName  The string of text / concept name to use.
#'@param schema        The database schema being used.
#'@param dbms          The target DBMS for SQL to be rendered in.
#'
#'@details Takes the aphroditeConceptName looks for synonyms and builds a list
#'of related concepts using the vocabulary hierarchies
#'
#'@return A list with two elements: a list of positive keywords found
#'(keywordlist_ALL), and a list of ignore keywords (ignorelist_ALL)
#'
#' @examples \dontrun{
#'
#'wordLists <- buildKeywordList(conn, aphrodite_concept_name, cdmSchema, dbms)
#'
#' }
#'
#'@export

getHbADisb <- function(){
	cohort <- read.table("diabetesCohort")
	pid <- unique(cohort$patient)
	pidToKeep <- rep(0,length(pid))
	for(i in 1:length(pid)){
	  dat <- cohort[cohort$patient==pid[i],]
	  drugFreq <- as.data.frame(table(as.character(dat$drugid)))
	  if(length(drugFreq$Freq)>1) next
	  pidToKeep[i] <- 1
	}
	remove(drugFreq,dat)
	pid <- pid[which(pidToKeep==1)]
	remove(pidToKeep)
	#SubCohort of patient who remained on drug with what they started 
	cohortSameStartEnd <- data.frame()
	for(i in 1:length(pid)){
	  dat <- cohort[cohort$patient==pid[i],]
	  pidStartDate <- as.Date(dat$pidindexdate[1])
	  pidEndDate <- as.Date(dat$pidindexdate[length(dat$patient)])
	  pidStartDrug <- as.character(dat$drugid[1])
	  pidEndDrug <- as.character(dat$drugid[1])
	  pidGender <- dat$gender[1]
	  pidAge <- dat$age[1]
	  dat2 <- cbind(pid[i],as.character(pidStartDate),as.character(pidEndDate),pidStartDrug,pidEndDrug,pidGender,pidAge)
	  cohortSameStartEnd <- rbind(cohortSameStartEnd,dat2)
	}
	colnames(cohortSameStartEnd) <- c("patient","pidindexdate","pidenddate","pidstartdrug","pidenddrug","pidgender","pidage")
	cohort <- cohort[!(cohort$patient %in% pid),]
	remove(pid,dat2,pidStartDate,pidEndDate,pidStartDrug,pidEndDrug,pidGender,pidAge,dat)

	#SubCase-1 These patients should move to other drug not before 30 days !
	pid <- unique(cohort$patient)
	cohortStartSwitch <- data.frame()
	for(i in 1:length(pid)){
	  dat <- cohort[cohort$patient==pid[i],]
	  pidStartDate <- as.Date(dat$pidindexdate[1])
	  pidEndDate <- as.Date(dat$pidindexdate[length(dat$patient)])
	  pidStartDrug <- as.character(dat$drugid[1])
	  pidGender <- dat$gender[1]
	  pidAge <- dat$age[1]
	  for(j in 2:length(dat$patient)){
		x <- as.numeric(grep(paste("^",dat$drugid[j-1],"$",sep=""),dat$drugid[j]))
		if(length(x)!=1) break
	  }
	  if(j==2){
		pidMovedDate <- as.Date(dat$pidindexdate[j])
		pidMovedDrug <- as.character(dat$drugid[j])
	  }else
	  {
		pidMovedDate <- as.Date(dat$pidindexdate[j])
		pidMovedDrug <- as.character(dat$drugid[j])
	  }
	  remove(j)
	  if(as.numeric(pidMovedDate-pidStartDate)<30) next #Change here if you wanat to cut loose this constraint
	  dat2 <- cbind(pid[i],as.character(pidStartDate),as.character(pidMovedDate),as.character(pidEndDate),pidStartDrug,pidMovedDrug,pidGender,pidAge)
	  cohortStartSwitch <- rbind(cohortStartSwitch,dat2)
	}
	colnames(cohortStartSwitch) <- c("patient","pidindexdate","pidmoveddate","pidenddate","pidstartdrug","pidmoveddrug","pidgender","pidage")
	remove(pid,dat2,pidStartDate,pidEndDate,pidStartDrug,pidMovedDrug,pidMovedDate,x,i,cohort,pidGender,pidAge,dat)
	# Dividing the patients into respective age groups ...
	cohortSameStartEnd$pidage <- as.numeric(as.character(cohortSameStartEnd$pidage))
	#ageGp-1 age <- 45
	cohortSameStartEndGpOne <- filter(cohortSameStartEnd,pidage<=45)
	# ageGp-2 46 < age < 75
	cohortSameStartEndGpTwo <- filter(cohortSameStartEnd,as.numeric(as.character(pidage))>=46&as.numeric(as.character(pidage))<=75)
	# ageGp-2 76 < age < infinity
	cohortSameStartEndGpThree <- filter(cohortSameStartEnd,as.numeric(as.character(pidage))>=76)

	cohortStartSwitchGpOne <- filter(cohortStartSwitch,as.numeric(as.character(pidage))<=45)
	cohortStartSwitchGpTwo <- filter(cohortStartSwitch,as.numeric(as.character(pidage))>=46&as.numeric(as.character(pidage))<=75)
	cohortStartSwitchGpThree <- filter(cohortStartSwitch,as.numeric(as.character(pidage))>=76)
	cohortStartSwitchGpOne$drugPath <- paste(cohortStartSwitchGpOne$pidstartdrug,cohortStartSwitchGpOne$pidmoveddrug,sep="->")
	cohortStartSwitchGpTwo$drugPath <- paste(cohortStartSwitchGpTwo$pidstartdrug,cohortStartSwitchGpTwo$pidmoveddrug,sep="->")
	cohortStartSwitchGpThree$drugPath <- paste(cohortStartSwitchGpThree$pidstartdrug,cohortStartSwitchGpThree$pidmoveddrug,sep="->")
	# Get the HbA1c for each of these groups ....
	measurements <- read.table("measurements")
	measurements$measureddate <- as.Date(measurements$measureddate)

	#------- FOR GROUP - 1 ------------------
	drug <- as.data.frame(table(cohortSameStartEndGpOne$pidstartdrug))
	drug <- drug[order(-drug$Freq),]
	drug <- filter(drug,Freq>=10)
	cohortSameStartEndDistGpOne <- data.frame()
	for(i in 1:length(drug$Var1)){
	  cohortSameStartEndDistGpOne <- rbind(cohortSameStartEndDistGpOne,cohortSameStartEndGpOne[cohortSameStartEndGpOne$pidstartdrug==drug$Var1[i],])
	}
	pid <- as.numeric(as.character(cohortSameStartEndDistGpOne$patient))
	pidLab <- data.frame()
	for(i in 1:length(pid)){
	  pidDat <- cohortSameStartEndGpOne[cohortSameStartEndGpOne$patient==pid[i],]
	  pidStart <- as.Date(pidDat$pidindexdate)
	  pidEnd <- as.Date(pidDat$pidenddate)
	  pidMes <- measurements[measurements$patient==pid[i],]
	  pidMes <- filter(pidMes,measureddate > as.Date(pidStart))
	  if(dim(pidMes)[1]==0) next
	  pidLab <- rbind(pidLab,cbind(rep(pid[i],length(pidMes$labval)),pidMes$labval))
	}
	colnames(pidLab) <- c("patient","labval")
	hbASameStartEndGpOne <- merge(pidLab,cohortSameStartEndDistGpOne,by="patient")
	p1 <- ggplot(hbASameStartEndGpOne, aes(x = pidstartdrug, y = labval, fill = pidstartdrug)) + geom_boxplot(outlier.colour = "grey") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1),text = element_text(size=14)) + labs(title="Patient Started With Drug and Remained Gp-1")

	#------- FOR GROUP - 2 ------------------
	drug <- as.data.frame(table(cohortSameStartEndGpTwo$pidstartdrug))
	drug <- drug[order(-drug$Freq),]
	drug <- filter(drug,Freq>=10)
	cohortSameStartEndDistGpTwo <- data.frame()
	for(i in 1:length(drug$Var1)){
	  cohortSameStartEndDistGpTwo <- rbind(cohortSameStartEndDistGpTwo,cohortSameStartEndGpTwo[cohortSameStartEndGpTwo$pidstartdrug==drug$Var1[i],])
	}
	pid <- as.numeric(as.character(cohortSameStartEndDistGpTwo$patient))
	pidLab <- data.frame()
	for(i in 1:length(pid)){
	  pidDat <- cohortSameStartEndGpTwo[cohortSameStartEndGpTwo$patient==pid[i],]
	  pidStart <- as.Date(pidDat$pidindexdate)
	  pidEnd <- as.Date(pidDat$pidenddate)
	  pidMes <- measurements[measurements$patient==pid[i],]
	  pidMes <- filter(pidMes,measureddate > as.Date(pidStart))
	  if(dim(pidMes)[1]==0) next
	  pidLab <- rbind(pidLab,cbind(rep(pid[i],length(pidMes$labval)),pidMes$labval))
	}
	colnames(pidLab) <- c("patient","labval")
	hbASameStartEndGpTwo <- merge(pidLab,cohortSameStartEndDistGpTwo,by="patient")
	p2 <- ggplot(hbASameStartEndGpTwo, aes(x = pidstartdrug, y = labval, fill = pidstartdrug)) + geom_boxplot()  + geom_boxplot(outlier.colour = "grey") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1),text = element_text(size=14)) + labs(title="Patient Started With Drug and Remained Gp-2")

	#------- FOR GROUP - 3 ------------------
	drug <- as.data.frame(table(cohortSameStartEndGpThree$pidstartdrug))
	drug <- drug[order(-drug$Freq),]
	drug <- filter(drug,Freq>=10)
	cohortSameStartEndDistGpThree <- data.frame()
	for(i in 1:length(drug$Var1)){
	  cohortSameStartEndDistGpThree <- rbind(cohortSameStartEndDistGpThree,cohortSameStartEndGpThree[cohortSameStartEndGpThree$pidstartdrug==drug$Var1[i],])
	}
	pid <- as.numeric(as.character(cohortSameStartEndDistGpThree$patient))
	pidLab <- data.frame()
	for(i in 1:length(pid)){
	  pidDat <- cohortSameStartEndGpThree[cohortSameStartEndGpThree$patient==pid[i],]
	  pidStart <- as.Date(pidDat$pidindexdate)
	  pidEnd <- as.Date(pidDat$pidenddate)
	  pidMes <- measurements[measurements$patient==pid[i],]
	  pidMes <- filter(pidMes,measureddate > as.Date(pidStart))
	  if(dim(pidMes)[1]==0) next
	  pidLab <- rbind(pidLab,cbind(rep(pid[i],length(pidMes$labval)),pidMes$labval))
	}
	colnames(pidLab) <- c("patient","labval")
	hbASameStartEndGpThree <- merge(pidLab,cohortSameStartEndDistGpThree,by="patient")
	p3 <- ggplot(hbASameStartEndGpThree, aes(x = pidstartdrug, y = labval, fill = pidstartdrug)) + geom_boxplot() + geom_boxplot(outlier.colour = "grey") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1),text = element_text(size=14)) + labs(title="Patient Started With Drug and Remained Gp-3")

	# For - Group - 1 Case Started with a drug and moved to other drug --------
	drug <- as.data.frame(table(cohortStartSwitchGpOne$drugPath))
	drug <- drug[order(-drug$Freq),]
	drug <- filter(drug,Freq>=10)
	cohortStartSwitchDistGpOne <- data.frame()
	for(i in 1:length(drug$Var1)){
	  cohortStartSwitchDistGpOne <- rbind(cohortStartSwitchDistGpOne,cohortStartSwitchGpOne[cohortStartSwitchGpOne$drugPath==drug$Var1[i],])
	}
	pid <- as.numeric(as.character(cohortStartSwitchDistGpOne$patient))
	pidLab <- data.frame()
	for(i in 1:length(pid)){
	  pidDat <- cohortStartSwitchGpOne[cohortStartSwitchGpOne$patient==pid[i],]
	  pidStart <- as.Date(pidDat$pidmoveddate)
	  pidMes <- measurements[measurements$patient==pid[i],]
	  pidMes <- filter(pidMes,measureddate > as.Date(pidStart))
	  if(dim(pidMes)[1]==0) next
	  pidLab <- rbind(pidLab,cbind(rep(pid[i],length(pidMes$labval)),pidMes$labval))
	}
	colnames(pidLab) <- c("patient","labval")
	hbAStartSwitchGpOne <- merge(pidLab,cohortStartSwitchDistGpOne,by="patient")
	p4 <- ggplot(hbAStartSwitchGpOne, aes(x = drugPath, y = labval, fill = drugPath)) + geom_boxplot() + geom_boxplot(outlier.colour = "grey") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1),text = element_text(size=14)) + labs(title="Patient Started With Drug and Moved Gp-1")

	# For - Group - 2 Case Started with a drug and moved to other --------
	drug <- as.data.frame(table(cohortStartSwitchGpTwo$drugPath))
	drug <- drug[order(-drug$Freq),]
	drug <- filter(drug,Freq>=10)
	cohortStartSwitchDistGpTwo <- data.frame()
	for(i in 1:length(drug$Var1)){
	  cohortStartSwitchDistGpTwo <- rbind(cohortStartSwitchDistGpTwo,cohortStartSwitchGpTwo[cohortStartSwitchGpTwo$drugPath==drug$Var1[i],])
	}
	pid <- as.numeric(as.character(cohortStartSwitchDistGpTwo$patient))
	pidLab <- data.frame()
	for(i in 1:length(pid)){
	  pidDat <- cohortStartSwitchGpTwo[cohortStartSwitchGpTwo$patient==pid[i],]
	  pidStart <- as.Date(pidDat$pidmoveddate)
	  pidMes <- measurements[measurements$patient==pid[i],]
	  pidMes <- filter(pidMes,measureddate > as.Date(pidStart))
	  if(dim(pidMes)[1]==0) next
	  pidLab <- rbind(pidLab,cbind(rep(pid[i],length(pidMes$labval)),pidMes$labval))
	}
	colnames(pidLab) <- c("patient","labval")
	hbAStartSwitchGpTwo <- merge(pidLab,cohortStartSwitchDistGpTwo,by="patient")
	p5 <- ggplot(hbAStartSwitchGpTwo, aes(x = drugPath, y = labval, fill = drugPath)) + geom_boxplot()+ theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1),text = element_text(size=14)) + labs(title="Patient Started With Drug and Moved Gp-2")

	# For - Group - 3 Case Started with a drug and moved to other drug --------
	drug <- as.data.frame(table(cohortStartSwitchGpThree$drugPath))
	drug <- drug[order(-drug$Freq),]
	drug <- filter(drug,Freq>=5) # In this case I didn't found more than 10 patient so changed it to 5. Make sure to change it to 10 in your analysis
	cohortStartSwitchDistGpThree <- data.frame()
	for(i in 1:length(drug$Var1)){
	  cohortStartSwitchDistGpThree <- rbind(cohortStartSwitchDistGpThree,cohortStartSwitchGpThree[cohortStartSwitchGpThree$drugPath==drug$Var1[i],])
	}
	pid <- as.numeric(as.character(cohortStartSwitchDistGpThree$patient))
	pidLab <- data.frame()
	for(i in 1:length(pid)){
	  pidDat <- cohortStartSwitchGpThree[cohortStartSwitchGpThree$patient==pid[i],]
	  pidStart <- as.Date(pidDat$pidmoveddate)
	  pidMes <- measurements[measurements$patient==pid[i],]
	  pidMes <- filter(pidMes,measureddate > as.Date(pidStart))
	  if(dim(pidMes)[1]==0) next
	  pidLab <- rbind(pidLab,cbind(rep(pid[i],length(pidMes$labval)),pidMes$labval))
	}
	colnames(pidLab) <- c("patient","labval")
	hbAStartSwitchGpThree <- merge(pidLab,cohortStartSwitchDistGpThree,by="patient")
	p6 <- ggplot(hbAStartSwitchGpThree, aes(x = drugPath, y = labval, fill = drugPath)) + geom_boxplot() + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1),text = element_text(size=14)) + labs(title="Patient Started With Drug and Moved Gp-3")

	# Ploting all together
	pdf("HbA1cPlots.pdf")
	print(p1)
	print(p2)
	print(p3)
	print(p4)
	print(p5)
	print(p6)
	dev.off()
	remove(cohortSameStartEnd,cohortSameStartEndDistGpOne,cohortSameStartEndDistGpTwo,cohortSameStartEndDistGpThree,cohortSameStartEndGpOne,cohortSameStartEndGpTwo,cohortSameStartEndGpThree,cohortStartSwitch,cohortStartSwitchDistGpOne,cohortStartSwitchDistGpTwo,cohortStartSwitchDistGpThree,cohortStartSwitchGpOne,cohortStartSwitchGpTwo,cohortStartSwitchGpThree,drug,hbASameStartEndGpOne,hbASameStartEndGpTwo,hbASameStartEndGpThree,hbAStartSwitchGpOne,hbAStartSwitchGpTwo,hbAStartSwitchGpThree,measurements,pidDat,pidLab,pidMes,i,p1,p2,p3,p4,p5,p6,pid,pidEnd,pidStart)
}

#'This function generates keyword and ignore lists based on the expansion of
#'concepts.
#'
#'@description Given any given concept_id or string of text this function
#'generates keyword and ignore lists based on the expansion of concepts (looking
#'at their synonyms).
#'
#'@param connection    The connection to the database server.
#'@param aphroditeConceptName  The string of text / concept name to use.
#'@param schema        The database schema being used.
#'@param dbms          The target DBMS for SQL to be rendered in.
#'
#'@details Takes the aphroditeConceptName looks for synonyms and builds a list
#'of related concepts using the vocabulary hierarchies
#'
#'@return A list with two elements: a list of positive keywords found
#'(keywordlist_ALL), and a list of ignore keywords (ignorelist_ALL)
#'
#' @examples \dontrun{
#'
#'wordLists <- buildKeywordList(conn, aphrodite_concept_name, cdmSchema, dbms)
#'
#' }
#'
#'@export

getSurvival <- function(){
#This will be called when sourcing aux_functions in the main script
#source('~/OHDSI-Final/OHDSI-FinalScripts/ggsurv.R')
#source('~/OHDSI-Final/OHDSI-FinalScripts/plotSurvival.R')
	cohort <- read.table("diabetesCohort")
	pid <- unique(cohort$patient)
	pidToKeep <- rep(0,length(pid))
	for(i in 1:length(pid)){
	  dat <- cohort[cohort$patient==pid[i],]
	  drugFreq <- as.data.frame(table(as.character(dat$drugid)))
	  if(length(drugFreq$Freq)>1) next
	  pidToKeep[i] <- 1
	  }
	remove(drugFreq,dat)
	pid <- pid[which(pidToKeep==1)]
	remove(pidToKeep)
	#SubCohort of patient who remained on drug with what they started 
	cohortSameStartEnd <- data.frame()
	for(i in 1:length(pid)){
	  dat <- cohort[cohort$patient==pid[i],]
	  pidStartDate <- as.Date(dat$pidindexdate[1])
	  pidEndDate <- as.Date(dat$pidindexdate[length(dat$patient)])
	  pidStartDrug <- as.character(dat$drugid[1])
	  pidEndDrug <- as.character(dat$drugid[1])
	  pidGender <- dat$gender[1]
	  pidAge <- dat$age[1]
	  dat2 <- cbind(pid[i],as.character(pidStartDate),as.character(pidEndDate),pidStartDrug,pidEndDrug,pidGender,pidAge)
	  cohortSameStartEnd <- rbind(cohortSameStartEnd,dat2)
	}
	colnames(cohortSameStartEnd) <- c("patient","pidindexdate","pidenddate","pidstartdrug","pidenddrug","pidgender","pidage")
	cohort <- cohort[!(cohort$patient %in% pid),]
	remove(pid,dat2,pidStartDate,pidEndDate,pidStartDrug,pidEndDrug,pidGender,pidAge,dat)

	#SubCase-1 These patients should move to other drug not before 30 days !
	pid <- unique(cohort$patient)
	cohortStartSwitch <- data.frame()
	for(i in 1:length(pid)){
	  dat <- cohort[cohort$patient==pid[i],]
	  pidStartDate <- as.Date(dat$pidindexdate[1])
	  pidEndDate <- as.Date(dat$pidindexdate[length(dat$patient)])
	  pidStartDrug <- as.character(dat$drugid[1])
	  pidGender <- dat$gender[1]
	  pidAge <- dat$age[1]
		for(j in 2:length(dat$patient)){
		  x <- as.numeric(grep(paste("^",dat$drugid[j-1],"$",sep=""),dat$drugid[j]))
		  if(length(x)!=1) break
		}
		if(j==2){
		  pidMovedDate <- as.Date(dat$pidindexdate[j])
		  pidMovedDrug <- as.character(dat$drugid[j])
		}else
		{
		  pidMovedDate <- as.Date(dat$pidindexdate[j])
		  pidMovedDrug <- as.character(dat$drugid[j])
		}
	  remove(j)
	  if(as.numeric(pidMovedDate-pidStartDate)<30) next #Change here if you want to cut loose the constraint
	  dat2 <- cbind(pid[i],as.character(pidStartDate),as.character(pidMovedDate),as.character(pidEndDate),pidStartDrug,pidMovedDrug,pidGender,pidAge)
	  cohortStartSwitch <- rbind(cohortStartSwitch,dat2)
	}
	colnames(cohortStartSwitch) <- c("patient","pidindexdate","pidmoveddate","pidenddate","pidstartdrug","pidmoveddrug","pidgender","pidage")
	remove(pid,dat2,pidStartDate,pidEndDate,pidStartDrug,pidMovedDrug,pidMovedDate,x,i,cohort,pidGender,pidAge,dat)

	# Dividing the patients into respective age groups ...
	# We need to divid these patient into three age group
	cohortSameStartEnd$pidage <- as.numeric(as.character(cohortSameStartEnd$pidage))
	#ageGp-1 age <- 45
	cohortSameStartEndGpOne <- filter(cohortSameStartEnd,pidage<=45)
	# ageGp-2 46 < age < 75
	cohortSameStartEndGpTwo <- filter(cohortSameStartEnd,as.numeric(as.character(pidage))>=46&as.numeric(as.character(pidage))<=75)
	# ageGp-2 76 < age < infinity
	cohortSameStartEndGpThree <- filter(cohortSameStartEnd,as.numeric(as.character(pidage))>=76)

	cohortStartSwitchGpOne <- filter(cohortStartSwitch,as.numeric(as.character(pidage))<=45)
	cohortStartSwitchGpTwo <- filter(cohortStartSwitch,as.numeric(as.character(pidage))>=46&as.numeric(as.character(pidage))<=75)
	cohortStartSwitchGpThree <- filter(cohortStartSwitch,as.numeric(as.character(pidage))>=76)
	cohortStartSwitchGpOne$drugPath <- paste(cohortStartSwitchGpOne$pidstartdrug,cohortStartSwitchGpOne$pidmoveddrug,sep="->")
	cohortStartSwitchGpTwo$drugPath <- paste(cohortStartSwitchGpTwo$pidstartdrug,cohortStartSwitchGpTwo$pidmoveddrug,sep="->")
	cohortStartSwitchGpThree$drugPath <- paste(cohortStartSwitchGpThree$pidstartdrug,cohortStartSwitchGpThree$pidmoveddrug,sep="->")

	# Get the HbA1c for each of these groups ....
	measurements <- read.table("measurements")
	measurements$measureddate <- as.Date(measurements$measureddate)

	#-------------- This is for group one patient only ---------
	survDatSameStartEndGpOne <- data.frame()
	pid <- as.numeric(as.character(cohortSameStartEndGpOne$patient))
	for(i in 1:length(pid)){
	  pidDat <- cohortSameStartEndGpOne[cohortSameStartEndGpOne$patient==pid[i],]
	  pidStart <- as.Date(pidDat$pidindexdate)
	  pidEnd <- as.Date(pidDat$pidenddate)
	  pidStartDrug <- as.character(pidDat$pidstartdrug)
	  pidEndDrug <- as.character(pidDat$pidenddrug)
	  pidAge <- as.numeric(as.character(pidDat$pidage))
	  pidGender <- as.numeric(as.character(pidDat$pidgender))
	  pidMes <- measurements[measurements$patient==pid[i],]
	  pidMes <- filter(pidMes,measureddate > as.Date(pidStart))
	  labRef <- as.numeric(c("7"))
	  if(dim(pidMes)[1]==0){
		timeToEvent <- as.numeric(as.Date(pidEnd)-as.Date(pidStart))
		event <- c("0")
	  }else
	  {
		pidMes <- pidMes[order(pidMes$measureddate),]
		pidEnd <- as.Date(pidMes$measureddate[length(pidMes$patient)])
		timeToEvent <- as.numeric(as.Date(pidEnd)-as.Date(pidStart))
		event <- c("1")
	  }
	  dat2 <- cbind(pid[i],pidStartDrug,pidEndDrug,timeToEvent,event,pidAge,pidGender)
	  survDatSameStartEndGpOne <- rbind(survDatSameStartEndGpOne,dat2)
	}                  
	colnames(survDatSameStartEndGpOne)[1] <- c("patient")
	#For group two patients ....
	survDatSameStartEndGpTwo <- data.frame()
	pid <- as.numeric(as.character(cohortSameStartEndGpTwo$patient))
	for(i in 1:length(pid)){
	  pidDat <- cohortSameStartEndGpTwo[cohortSameStartEndGpTwo$patient==pid[i],]
	  pidStart <- as.Date(pidDat$pidindexdate)
	  pidEnd <- as.Date(pidDat$pidenddate)
	  pidStartDrug <- as.character(pidDat$pidstartdrug)
	  pidEndDrug <- as.character(pidDat$pidenddrug)
	  pidAge <- as.numeric(as.character(pidDat$pidage))
	  pidGender <- as.numeric(as.character(pidDat$pidgender))
	  pidMes <- measurements[measurements$patient==pid[i],]
	  pidMes <- filter(pidMes,measureddate > as.Date(pidStart))
	  labRef <- as.numeric(c("7.5"))
	  if(dim(pidMes)[1]==0){
		timeToEvent <- as.numeric(as.Date(pidEnd)-as.Date(pidStart))
		event <- c("0")
	  }else
	  {
		pidMes <- pidMes[order(pidMes$measureddate),]
		pidEnd <- as.Date(pidMes$measureddate[length(pidMes$patient)])
		timeToEvent <- as.numeric(as.Date(pidEnd)-as.Date(pidStart))
		event <- c("1")
	  }
	  dat2 <- cbind(pid[i],pidStartDrug,pidEndDrug,timeToEvent,event,pidAge,pidGender)
	  survDatSameStartEndGpTwo <- rbind(survDatSameStartEndGpTwo,dat2)
	} 
	colnames(survDatSameStartEndGpTwo)[1] <- c("patient")
	#For GroupThree patient
	survDatSameStartEndGpThree <- data.frame()
	pid <- as.numeric(as.character(cohortSameStartEndGpThree$patient))
	for(i in 1:length(pid)){
	  pidDat <- cohortSameStartEndGpThree[cohortSameStartEndGpThree$patient==pid[i],]
	  pidStart <- as.Date(pidDat$pidindexdate)
	  pidEnd <- as.Date(pidDat$pidenddate)
	  pidStartDrug <- as.character(pidDat$pidstartdrug)
	  pidEndDrug <- as.character(pidDat$pidenddrug)
	  pidAge <- as.numeric(as.character(pidDat$pidage))
	  pidGender <- as.numeric(as.character(pidDat$pidgender))
	  pidMes <- measurements[measurements$patient==pid[i],]
	  pidMes <- filter(pidMes,measureddate > as.Date(pidStart))
	  labRef <- as.numeric(c("8"))
	  if(dim(pidMes)[1]==0){
		timeToEvent <- as.numeric(as.Date(pidEnd)-as.Date(pidStart))
		event <- c("0")
	  }else
	  {
		pidMes <- pidMes[order(pidMes$measureddate),]
		pidEnd <- as.Date(pidMes$measureddate[length(pidMes$patient)])
		timeToEvent <- as.numeric(as.Date(pidEnd)-as.Date(pidStart))
		event <- c("1")
	  }
	  dat2 <- cbind(pid[i],pidStartDrug,pidEndDrug,timeToEvent,event,pidAge,pidGender)
	  survDatSameStartEndGpThree <- rbind(survDatSameStartEndGpThree,dat2)
	}
	colnames(survDatSameStartEndGpThree)[1] <- c("patient")

	#------------ For the guys who moved to other drugs / changed drug -------------
	#For GroupOne patient
	survDatStartSwitchGpOne <- data.frame()
	pid <- as.numeric(as.character(cohortStartSwitchGpOne$patient))
	for(i in 1:length(pid)){
	  pidDat <- cohortStartSwitchGpOne[cohortStartSwitchGpOne$patient==pid[i],]
	  pidStart <- as.Date(pidDat$pidindexdate)
	  pidMoved <- as.Date(pidDat$pidmoveddate)
	  pidEnd <- as.Date(pidDat$pidenddate)
	  pidAge <- as.numeric(as.character(pidDat$pidage))
	  pidGender <- as.numeric(as.character(pidDat$pidgender))
	  pidStartDrug <- as.character(pidDat$pidstartdrug)
	  pidEndDrug <- as.character(pidDat$pidmoveddrug)
	  pidMes <- measurements[measurements$patient==pid[i],]
	  pidMes <- filter(pidMes,measureddate > as.Date(pidMoved))
	  labRef <- as.numeric(c("7"))
	  if(dim(pidMes)[1]==0){
		timeToEvent <- as.numeric(as.Date(pidEnd)-as.Date(pidStart))
		event <- c("0")
	  }else
	  {
		pidMes <- pidMes[order(pidMes$measureddate),]
		pidEnd <- as.Date(pidMes$measureddate[length(pidMes$patient)])
		timeToEvent <- as.numeric(as.Date(pidEnd)-as.Date(pidStart))
		event <- c("1")
	  }
	  dat2 <- cbind(pid[i],pidStartDrug,pidEndDrug,timeToEvent,event,pidAge,pidGender)
	  survDatStartSwitchGpOne <- rbind(survDatStartSwitchGpOne,dat2)
	}
	colnames(survDatStartSwitchGpOne)[1] <- c("patient")
	#For GroupTwo patient
	survDatStartSwitchGpTwo <- data.frame()
	pid <- as.numeric(as.character(cohortStartSwitchGpTwo$patient))
	for(i in 1:length(pid)){
	  pidDat <- cohortStartSwitchGpTwo[cohortStartSwitchGpTwo$patient==pid[i],]
	  pidStart <- as.Date(pidDat$pidindexdate)
	  pidMoved <- as.Date(pidDat$pidmoveddate)
	  pidEnd <- as.Date(pidDat$pidenddate)
	  pidAge <- as.numeric(as.character(pidDat$pidage))
	  pidGender <- as.numeric(as.character(pidDat$pidgender))
	  pidStartDrug <- as.character(pidDat$pidstartdrug)
	  pidEndDrug <- as.character(pidDat$pidmoveddrug)
	  pidMes <- measurements[measurements$patient==pid[i],]
	  pidMes <- filter(pidMes,measureddate > as.Date(pidMoved))
	  labRef <- as.numeric(c("7.5"))
	  if(dim(pidMes)[1]==0){
		timeToEvent <- as.numeric(as.Date(pidEnd)-as.Date(pidStart))
		event <- c("0")
	  }else
	  {
		pidMes <- pidMes[order(pidMes$measureddate),]
		pidEnd <- as.Date(pidMes$measureddate[length(pidMes$patient)])
		timeToEvent <- as.numeric(as.Date(pidEnd)-as.Date(pidStart))
		event <- c("1")
	  }
	  dat2 <- cbind(pid[i],pidStartDrug,pidEndDrug,timeToEvent,event,pidAge,pidGender)
	  survDatStartSwitchGpTwo <- rbind(survDatStartSwitchGpTwo,dat2)
	}
	colnames(survDatStartSwitchGpTwo)[1] <- c("patient")
	#For GroupThree patient
	survDatStartSwitchGpThree <- data.frame()
	pid <- as.numeric(as.character(cohortStartSwitchGpThree$patient))
	for(i in 1:length(pid)){
	  pidDat <- cohortStartSwitchGpThree[cohortStartSwitchGpThree$patient==pid[i],]
	  pidStart <- as.Date(pidDat$pidindexdate)
	  pidMoved <- as.Date(pidDat$pidmoveddate)
	  pidEnd <- as.Date(pidDat$pidenddate)
	  pidAge <- as.numeric(as.character(pidDat$pidage))
	  pidGender <- as.numeric(as.character(pidDat$pidgender))
	  pidStartDrug <- as.character(pidDat$pidstartdrug)
	  pidEndDrug <- as.character(pidDat$pidmoveddrug)
	  pidMes <- measurements[measurements$patient==pid[i],]
	  pidMes <- filter(pidMes,measureddate > as.Date(pidMoved))
	  labRef <- as.numeric(c("8"))
	  if(dim(pidMes)[1]==0){
		timeToEvent <- as.numeric(as.Date(pidEnd)-as.Date(pidStart))
		event <- c("0")
	  }else
	  {
		pidMes <- pidMes[order(pidMes$measureddate),]
		pidEnd <- as.Date(pidMes$measureddate[length(pidMes$patient)])
		timeToEvent <- as.numeric(as.Date(pidEnd)-as.Date(pidStart))
		event <- c("1")
	  }
	  dat2 <- cbind(pid[i],pidStartDrug,pidEndDrug,timeToEvent,event,pidAge,pidGender)
	  survDatStartSwitchGpThree <- rbind(survDatStartSwitchGpThree,dat2)
	}
	colnames(survDatStartSwitchGpThree)[1] <- c("patient")
	#Get the final aggregated data for each group
	#Group-1
	survDatSameStartEndGpOne$drugPath <- paste(survDatSameStartEndGpOne$pidStartDrug,survDatSameStartEndGpOne$pidEndDrug,sep="->")
	survDatStartSwitchGpOne$drugPath <- paste(survDatStartSwitchGpOne$pidStartDrug,survDatStartSwitchGpOne$pidEndDrug,sep="->")
	survDatGpOne <- rbind(survDatSameStartEndGpOne,survDatStartSwitchGpOne)
	#Group-2
	survDatSameStartEndGpTwo$drugPath <- paste(survDatSameStartEndGpTwo$pidStartDrug,survDatSameStartEndGpTwo$pidEndDrug,sep="->")
	survDatStartSwitchGpTwo$drugPath <- paste(survDatStartSwitchGpTwo$pidStartDrug,survDatStartSwitchGpTwo$pidEndDrug,sep="->")
	survDatGpTwo <- rbind(survDatSameStartEndGpTwo,survDatStartSwitchGpTwo)
	#Group-3
	survDatSameStartEndGpThree$drugPath <- paste(survDatSameStartEndGpThree$pidStartDrug,survDatSameStartEndGpThree$pidEndDrug,sep="->")
	survDatStartSwitchGpThree$drugPath <- paste(survDatStartSwitchGpThree$pidStartDrug,survDatStartSwitchGpThree$pidEndDrug,sep="->")
	survDatGpThree <- rbind(survDatSameStartEndGpThree,survDatStartSwitchGpThree)
	remove(cohortSameStartEnd,cohortSameStartEndGpOne,cohortSameStartEndGpTwo,cohortSameStartEndGpThree,cohortStartSwitch,cohortStartSwitchGpOne,cohortStartSwitchGpTwo,cohortStartSwitchGpThree,dat2,measurements,pidDat,pidMes,survDatStartSwitchGpThree,survDatStartSwitchGpTwo,survDatStartSwitchGpOne,survDatSameStartEndGpOne,survDatSameStartEndGpTwo,survDatSameStartEndGpThree,event,i,labRef,pid,pidAge,pidEnd,pidEndDrug,pidGender,pidMoved,pidStart,pidStartDrug,timeToEvent)

	#----------- Survival and plotting -------------
	#HR plot - logrank nonparametric
	#Group-1
	pidDrug <- as.data.frame(table(survDatGpOne$drugPath))
	pidDrug <- pidDrug[order(-pidDrug$Freq),]
	pidDrug <- pidDrug[which(pidDrug$Freq>=15),]
	colnames(pidDrug) <- c("drugComb","Freq")
	survDat <- list()
	for(i in 1:length(pidDrug$drugComb)){
	  survDat[[i]] <- survDatGpOne[survDatGpOne$drugPath==pidDrug$drugComb[i],]
	}
	#pairwise survival
	pair <- t(combn(length(survDat),2))
	survDatPair <- list()
	for(i in 1:length(pair[,1])){
	  survDatPair[[i]] <- rbind(survDat[[pair[i,1]]],survDat[[pair[i,2]]])
	}
	#Survival plots
	survFit <- list()
	for(i in 1:length(survDatPair)){
	  dat <- survDatPair[[i]]
	  dat$timeToEvent <- as.numeric(as.character(dat$timeToEvent))
	  dat$event <- as.numeric(as.character(dat$event))
	  survFit[[i]] <- survfit(Surv(timeToEvent,event)~drugPath,data=dat)
	}
	#Survival diff - log rank non parametric ...
	survDiff <- list()
	for(i in 1:length(survDatPair)){
	  dat <- survDatPair[[i]]
	  dat$timeToEvent <- as.numeric(as.character(dat$timeToEvent))
	  dat$event <- as.numeric(as.character(dat$event))
	  survDiff[[i]] <- survdiff(Surv(timeToEvent,event)~drugPath,data=dat)
	}
	#Log-Rank test un-adjusted
	survLogRank <- data.frame()
	for(i in 1:length(survDatPair)){
	  dat22 <- survDatPair[[i]]
	  dat22$timeToEvent <- as.numeric(as.character(dat22$timeToEvent))
	  dat22$event <- as.numeric(as.character(dat22$event))
	  fit <- survdiff(Surv(timeToEvent,event)~drugPath,data=dat22)
	  pVal = 1 - pchisq(fit$chisq, length(fit$n) - 1)
	  HR = (fit$obs[2]/fit$exp[2])/(fit$obs[1]/fit$exp[1])
	  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/fit$exp[2]+1/fit$exp[1]))
	  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/fit$exp[2]+1/fit$exp[1]))
	  dat <- cbind(names(fit$n[1]),names(fit$n[2]),pVal,HR,up95,low95)
	  survLogRank <- rbind(survLogRank,dat)
	}
	colnames(survLogRank)[1] <- c("drugPath1Ref")
	colnames(survLogRank)[2] <- c("drugPath2")
	survLogRank$pVal <- format(round(as.numeric(as.character(survLogRank$pVal)), 4), nsmall = 4)
	survLogRank$HR <- format(round(as.numeric(as.character(survLogRank$HR)), 2), nsmall = 2)
	survLogRank$up95 <- format(round(as.numeric(as.character(survLogRank$up95)), 2), nsmall = 2)
	survLogRank$low95 <- format(round(as.numeric(as.character(survLogRank$low95)), 2), nsmall = 2)
	label <- as.character((1:length(survLogRank$drugPath1Ref)))
	### Only For plotting ####
	dat <- cbind(label,survLogRank[,c("HR","low95","up95")])
	dat$label <- factor(dat$label, levels=dat$label)
	dat$HR <- as.numeric(dat$HR)
	dat$up95 <- as.numeric(dat$up95)
	dat$low95 <- as.numeric(dat$low95)
	p1 <- ggplot(dat, aes(x=label, y=HR, ymin=low95, ymax=up95)) + geom_pointrange() + coord_flip() + geom_point(color="maroon",shape=15,size=5) + geom_hline(yintercept=1,linetype="dashed")+theme_bw()+ theme(text = element_text(size=14)) + labs(title="HR - Patient Started with a Drug and Remained Gp-1")
	#printing survlogRank
	colnames(survLogRank)[1] <- c("RefPath")
	colnames(survLogRank)[2] <- c("Path")
	survLogRank$RefPath <- gsub("drugPath=","",survLogRank$RefPath)
	survLogRank$Path <- gsub("drugPath=","",survLogRank$Path)
	survLogRank <- survLogRank[rev(rownames(survLogRank)),]
	p1Table <- survLogRank
	#grid.table(survLogRank)
	#dev.off()
	#########################  
	remove(dat,dat22,pair,pidDrug,fit,HR,low95,up95,pVal,label,i)

	#Group-2
	pidDrugGpTwo <- as.data.frame(table(survDatGpTwo$drugPath))
	pidDrugGpTwo <- pidDrugGpTwo[order(-pidDrugGpTwo$Freq),]
	pidDrugGpTwo <- pidDrugGpTwo[which(pidDrugGpTwo$Freq>=25),]
	colnames(pidDrugGpTwo) <- c("drugComb","Freq")
	survDatPlotGpTwo <- list()
	for(i in 1:length(pidDrugGpTwo$drugComb)){
	  survDatPlotGpTwo[[i]] <- survDatGpTwo[survDatGpTwo$drugPath==pidDrugGpTwo$drugComb[i],]
	}
	#pairwise survival
	pairGpTwo <- t(combn(length(survDatPlotGpTwo),2))
	survDatPairGpTwo <- list()
	for(i in 1:length(pairGpTwo[,1])){
	  survDatPairGpTwo[[i]] <- rbind(survDatPlotGpTwo[[pairGpTwo[i,1]]],survDatPlotGpTwo[[pairGpTwo[i,2]]])
	}
	#Survival plots
	survFitGpTwo <- list()
	for(i in 1:length(survDatPairGpTwo)){
	  dat <- survDatPairGpTwo[[i]]
	  dat$timeToEvent <- as.numeric(as.character(dat$timeToEvent))
	  dat$event <- as.numeric(as.character(dat$event))
	  survFitGpTwo[[i]] <- survfit(Surv(timeToEvent,event)~drugPath,data=dat)
	}
	#Survival diff - log rank non parametric ...
	survDiffGpTwo <- list()
	for(i in 1:length(survDatPairGpTwo)){
	  dat <- survDatPairGpTwo[[i]]
	  dat$timeToEvent <- as.numeric(as.character(dat$timeToEvent))
	  dat$event <- as.numeric(as.character(dat$event))
	  survDiffGpTwo[[i]] <- survdiff(Surv(timeToEvent,event)~drugPath,data=dat)
	}
	#Log-Rank test un-adjusted
	survLogRankGpTwo <- data.frame()
	for(i in 1:length(survDatPairGpTwo)){
	  dat22 <- survDatPairGpTwo[[i]]
	  dat22$timeToEvent <- as.numeric(as.character(dat22$timeToEvent))
	  dat22$event <- as.numeric(as.character(dat22$event))
	  fit <- survdiff(Surv(timeToEvent,event)~drugPath,data=dat22)
	  pVal = 1 - pchisq(fit$chisq, length(fit$n) - 1)
	  HR = (fit$obs[2]/fit$exp[2])/(fit$obs[1]/fit$exp[1])
	  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/fit$exp[2]+1/fit$exp[1]))
	  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/fit$exp[2]+1/fit$exp[1]))
	  dat <- cbind(names(fit$n[1]),names(fit$n[2]),pVal,HR,up95,low95)
	  survLogRankGpTwo <- rbind(survLogRankGpTwo,dat)
	}
	colnames(survLogRankGpTwo)[1] <- c("drugPath1Ref")
	colnames(survLogRankGpTwo)[2] <- c("drugPath2")
	survLogRankGpTwo$pVal <- format(round(as.numeric(as.character(survLogRankGpTwo$pVal)), 4), nsmall = 4)
	survLogRankGpTwo$HR <- format(round(as.numeric(as.character(survLogRankGpTwo$HR)), 2), nsmall = 2)
	survLogRankGpTwo$up95 <- format(round(as.numeric(as.character(survLogRankGpTwo$up95)), 2), nsmall = 2)
	survLogRankGpTwo$low95 <- format(round(as.numeric(as.character(survLogRankGpTwo$low95)), 2), nsmall = 2)
	label <- as.character((1:length(survLogRankGpTwo$drugPath1Ref)))
	### Only For plotting ####
	dat <- cbind(label,survLogRankGpTwo[,c("HR","low95","up95")])
	dat$label <- factor(dat$label, levels=dat$label)
	dat$HR <- as.numeric(dat$HR)
	dat$up95 <- as.numeric(dat$up95)
	dat$low95 <- as.numeric(dat$low95)
	#dat$age <- as.numeric(as.numeric((format(as.Date(dat$pidStartIndexDate),format="%Y")))-as.numeric(as.character(dat$pidBirthYear)))
	p2 <- ggplot(dat, aes(x=label, y=HR, ymin=low95, ymax=up95)) + geom_pointrange() + coord_flip() + geom_point(color="maroon",shape=15,size=5) + geom_hline(yintercept=1,linetype="dashed")+theme_bw()+theme(text = element_text(size=14))
	colnames(survLogRankGpTwo)[1] <- c("RefPath")
	colnames(survLogRankGpTwo)[2] <- c("Path")
	survLogRankGpTwo$RefPath <- gsub("drugPath=","",survLogRankGpTwo$RefPath)
	survLogRankGpTwo$Path <- gsub("drugPath=","",survLogRankGpTwo$Path)
	survLogRankGpTwo <- survLogRankGpTwo[rev(rownames(survLogRankGpTwo)),]
	p2Table <- survLogRankGpTwo
	#pdf("survLogRankGpTwo.pdf", height=length(survLogRankGpTwo$RefPath), width=8)
	#grid.table(survLogRankGpTwo)
	#dev.off()
	#########################  

	#Group-3
	pidDrugGpThree <- as.data.frame(table(survDatGpThree$drugPath))
	pidDrugGpThree <- pidDrugGpThree[order(-pidDrugGpThree$Freq),]
	pidDrugGpThree <- pidDrugGpThree[which(pidDrugGpThree$Freq>=10),]
	colnames(pidDrugGpThree) <- c("drugComb","Freq")
	survDatPlotGpThree <- list()
	for(i in 1:length(pidDrugGpThree$drugComb)){
	  survDatPlotGpThree[[i]] <- survDatGpThree[survDatGpThree$drugPath==pidDrugGpThree$drugComb[i],]
	}
	#pairwise survival
	pairGpThree <- t(combn(length(survDatPlotGpThree),2))
	survDatPairGpThree <- list()
	for(i in 1:length(pairGpThree[,1])){
	  survDatPairGpThree[[i]] <- rbind(survDatPlotGpThree[[pairGpThree[i,1]]],survDatPlotGpThree[[pairGpThree[i,2]]])
	}
	#Survival plots
	survFitGpThree <- list()
	for(i in 1:length(survDatPairGpThree)){
	  dat <- survDatPairGpThree[[i]]
	  dat$timeToEvent <- as.numeric(as.character(dat$timeToEvent))
	  dat$event <- as.numeric(as.character(dat$event))
	  survFitGpThree[[i]] <- survfit(Surv(timeToEvent,event)~drugPath,data=dat)
	}
	#Survival diff - log rank non parametric ...
	survDiffGpThree <- list()
	for(i in 1:length(survDatPairGpThree)){
	  dat <- survDatPairGpThree[[i]]
	  dat$timeToEvent <- as.numeric(as.character(dat$timeToEvent))
	  dat$event <- as.numeric(as.character(dat$event))
	  survDiffGpThree[[i]] <- survdiff(Surv(timeToEvent,event)~drugPath,data=dat)
	}
	#Log-Rank test un-adjusted
	survLogRankGpThree <- data.frame()
	for(i in 1:length(survDatPairGpThree)){
	  dat22 <- survDatPairGpThree[[i]]
	  dat22$timeToEvent <- as.numeric(as.character(dat22$timeToEvent))
	  dat22$event <- as.numeric(as.character(dat22$event))
	  fit <- survdiff(Surv(timeToEvent,event)~drugPath,data=dat22)
	  pVal = 1 - pchisq(fit$chisq, length(fit$n) - 1)
	  HR = (fit$obs[2]/fit$exp[2])/(fit$obs[1]/fit$exp[1])
	  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/fit$exp[2]+1/fit$exp[1]))
	  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/fit$exp[2]+1/fit$exp[1]))
	  dat <- cbind(names(fit$n[1]),names(fit$n[2]),pVal,HR,up95,low95)
	  survLogRankGpThree <- rbind(survLogRankGpThree,dat)
	}
	colnames(survLogRankGpThree)[1] <- c("drugPath1Ref")
	colnames(survLogRankGpThree)[2] <- c("drugPath2")
	survLogRankGpThree$pVal <- format(round(as.numeric(as.character(survLogRankGpThree$pVal)), 4), nsmall = 4)
	survLogRankGpThree$HR <- format(round(as.numeric(as.character(survLogRankGpThree$HR)), 2), nsmall = 2)
	survLogRankGpThree$up95 <- format(round(as.numeric(as.character(survLogRankGpThree$up95)), 2), nsmall = 2)
	survLogRankGpThree$low95 <- format(round(as.numeric(as.character(survLogRankGpThree$low95)), 2), nsmall = 2)
	label <- as.character((1:length(survLogRankGpThree$drugPath1Ref)))
	### Only For plotting ####
	dat <- cbind(label,survLogRankGpThree[,c("HR","low95","up95")])
	dat$label <- factor(dat$label, levels=dat$label)
	dat$HR <- as.numeric(dat$HR)
	dat$up95 <- as.numeric(dat$up95)
	dat$low95 <- as.numeric(dat$low95)
	#dat$age <- as.numeric(as.numeric((format(as.Date(dat$pidStartIndexDate),format="%Y")))-as.numeric(as.character(dat$pidBirthYear)))
	p3 <- ggplot(dat, aes(x=label, y=HR, ymin=low95, ymax=up95)) + geom_pointrange() + coord_flip() + geom_point(color="maroon",shape=15,size=5) + geom_hline(yintercept=1,linetype="dashed")+theme_bw()+theme(text = element_text(size=14))
	colnames(survLogRankGpThree)[1] <- c("RefPath")
	colnames(survLogRankGpThree)[2] <- c("Path")
	survLogRankGpThree$RefPath <- gsub("drugPath=","",survLogRankGpThree$RefPath)
	survLogRankGpThree$Path <- gsub("drugPath=","",survLogRankGpThree$Path)
	survLogRankGpThree <- survLogRankGpThree[rev(rownames(survLogRankGpThree)),]
	p3Table <- survLogRankGpThree
	#########################  
	pdf("survivalPlots.pdf")
	print(p1)
	#grid.table(p1Table)
	for(i in 1:length(survFit)){
	  plotSurvival(survFit[[i]],survDiff[[i]],"Group-1")
	}
	print(p2)
	for(i in 1:length(survFitGpTwo)){ # take it into the pdf loop
	  plotSurvival(survFitGpTwo[[i]],survDiffGpTwo[[i]],"Group-2")
	}
	print(p3)
	for(i in 1:length(survFitGpThree)){
	  plotSurvival(survFitGpThree[[i]],survDiffGpThree[[i]],"Group-3")
	}
	dev.off()

	tables <- list(p1Table,p2Table,p3Table)
	for(i in 1:length(tables)){
	  pdf(paste("survivalTableGp",i,".pdf",sep=""),length(tables[[i]]$RefPath),width=10)
	  grid.table(tables[[i]])
	  dev.off()
	}
	remove(dat,dat22,p1Table,p2Table,p3Table,pairGpTwo,pairGpThree,pidDrugGpThree,pidDrugGpTwo,survDatGpOne,survDatGpThree,survLogRank,survLogRankGpThree,survLogRankGpTwo,fit,HR,i,label,low95,p1,p2,p3,pVal,survDat,survDatPair,survDatPairGpThree,survDatPairGpTwo,survDatPlotGpThree,survDatPlotGpTwo,survDiff,survDiffGpThree,survDiffGpTwo,survFit,survFitGpThree,survFitGpTwo,tables,up95,survDatGpTwo)
}
#'This function generates keyword and ignore lists based on the expansion of
#'concepts.
#'
#'@description Given any given concept_id or string of text this function
#'generates keyword and ignore lists based on the expansion of concepts (looking
#'at their synonyms).
#'
#'@param connection    The connection to the database server.
#'@param aphroditeConceptName  The string of text / concept name to use.
#'@param schema        The database schema being used.
#'@param dbms          The target DBMS for SQL to be rendered in.
#'
#'@details Takes the aphroditeConceptName looks for synonyms and builds a list
#'of related concepts using the vocabulary hierarchies
#'
#'@return A list with two elements: a list of positive keywords found
#'(keywordlist_ALL), and a list of ignore keywords (ignorelist_ALL)
#'
#' @examples \dontrun{
#'
#'wordLists <- buildKeywordList(conn, aphrodite_concept_name, cdmSchema, dbms)
#'
#' }
#'
#'@export
etTreatmentPath <- function(){
	  cohort <- read.table("diabetesCohort")
	  pid <- unique(cohort$patient)
	  pidToKeep <- rep(0,length(pid))
	  for(i in 1:length(pid)){
		dat <- cohort[cohort$patient==pid[i],]
		drugFreq <- as.data.frame(table(as.character(dat$drugid)))
		if(length(drugFreq$Freq)>1) next
		pidToKeep[i] <- 1
	  }
	  remove(drugFreq,dat)
	  pid <- pid[which(pidToKeep==1)]
	  remove(pidToKeep)
	  #SubCohort of patient who remained on drug with what they started 
	  cohortSameStartEnd <- data.frame()
	  for(i in 1:length(pid)){
		dat <- cohort[cohort$patient==pid[i],]
		pidStartDate <- as.Date(dat$pidindexdate[1])
		pidEndDate <- as.Date(dat$pidindexdate[length(dat$patient)])
		pidStartDrug <- as.character(dat$drugid[1])
		pidEndDrug <- as.character(dat$drugid[1])
		pidGender <- dat$gender[1]
		pidAge <- dat$age[1]
		dat2 <- cbind(pid[i],as.character(pidStartDate),as.character(pidEndDate),pidStartDrug,pidEndDrug,pidGender,pidAge)
		cohortSameStartEnd <- rbind(cohortSameStartEnd,dat2)
	  }
	  colnames(cohortSameStartEnd) <- c("patient","pidindexdate","pidenddate","pidstartdrug","pidenddrug","pidgender","pidage")
	  #Case-2 Get the patient who started with one drug and move to second within at least 30 days ...
	  cohort <- cohort[!(cohort$patient %in% pid),] #these are the patient who moved to other drugs
	  remove(pid,dat2,pidStartDate,pidEndDate,pidStartDrug,pidEndDrug,pidGender,pidAge,dat)
	  #SubCase-1 These patients should move to other drug not before 30 days !
	  pid <- unique(cohort$patient)
	  cohortStartSwitch <- data.frame()
	  for(i in 1:length(pid)){
		dat <- cohort[cohort$patient==pid[i],]
		pidStartDate <- as.Date(dat$pidindexdate[1])
		pidEndDate <- as.Date(dat$pidindexdate[length(dat$patient)])
		pidStartDrug <- as.character(dat$drugid[1])
		pidGender <- dat$gender[1]
		pidAge <- dat$age[1]
		for(j in 2:length(dat$patient)){
		  x <- as.numeric(grep(paste("^",dat$drugid[j-1],"$",sep=""),dat$drugid[j]))
		  if(length(x)!=1) break
		}
		if(j==2){
		  pidMovedDate <- as.Date(dat$pidindexdate[j])
		  pidMovedDrug <- as.character(dat$drugid[j])
		}else
		{
		  pidMovedDate <- as.Date(dat$pidindexdate[j])
		  pidMovedDrug <- as.character(dat$drugid[j])
		}
		remove(j)
		if(as.numeric(pidMovedDate-pidStartDate)<30) next #Change here if you wanat to cut loose the constrain
		dat2 <- cbind(pid[i],as.character(pidStartDate),as.character(pidMovedDate),as.character(pidEndDate),pidStartDrug,pidMovedDrug,pidGender,pidAge)
		cohortStartSwitch <- rbind(cohortStartSwitch,dat2)
	  }
	  colnames(cohortStartSwitch) <- c("patient","pidindexdate","pidmoveddate","pidenddate","pidstartdrug","pidmoveddrug","pidgender","pidage")
	  remove(pid,dat2,pidStartDate,pidEndDate,pidStartDrug,pidMovedDrug,pidMovedDate,x,i,cohort,pidGender,pidAge,dat)
	  # Dividing the patients into respective age groups ...
	  cohortSameStartEnd$pidage <- as.numeric(as.character(cohortSameStartEnd$pidage))
	  cohortSameStartEndGpOne <- filter(cohortSameStartEnd,pidage<=45)
	  cohortSameStartEndGpOne$drugPath <- paste(cohortSameStartEndGpOne$pidstartdrug,cohortSameStartEndGpOne$pidenddrug,sep="->")
	  cohortSameStartEndGpTwo <- filter(cohortSameStartEnd,as.numeric(as.character(pidage))>=46&as.numeric(as.character(pidage))<=75)
	  cohortSameStartEndGpTwo$drugPath <- paste(cohortSameStartEndGpTwo$pidstartdrug,cohortSameStartEndGpTwo$pidenddrug,sep="->")
	  cohortSameStartEndGpThree <- filter(cohortSameStartEnd,as.numeric(as.character(pidage))>=76)
	  cohortSameStartEndGpThree$drugPath <- paste(cohortSameStartEndGpThree$pidstartdrug,cohortSameStartEndGpThree$pidenddrug,sep="->")
	  cohortStartSwitchGpOne <- filter(cohortStartSwitch,as.numeric(as.character(pidage))<=45)
	  cohortStartSwitchGpTwo <- filter(cohortStartSwitch,as.numeric(as.character(pidage))>=46&as.numeric(as.character(pidage))<=75)
	  cohortStartSwitchGpThree <- filter(cohortStartSwitch,as.numeric(as.character(pidage))>=76)
	  cohortStartSwitchGpOne$drugPath <- paste(cohortStartSwitchGpOne$pidstartdrug,cohortStartSwitchGpOne$pidmoveddrug,sep="->")
	  cohortStartSwitchGpTwo$drugPath <- paste(cohortStartSwitchGpTwo$pidstartdrug,cohortStartSwitchGpTwo$pidmoveddrug,sep="->")
	  cohortStartSwitchGpThree$drugPath <- paste(cohortStartSwitchGpThree$pidstartdrug,cohortStartSwitchGpThree$pidmoveddrug,sep="->")
	  #Group-1
	  drugsSameStartEndGpOne <- as.data.frame(table(cohortSameStartEndGpOne$drugPath))
	  drugsSameStartEndGpOne <- drugsSameStartEndGpOne[order(-drugsSameStartEndGpOne$Freq),]
	  dat <- data.frame(do.call('rbind', strsplit(as.character(drugsSameStartEndGpOne$Var1),'->',fixed=TRUE)))
	  dat <- cbind(dat,drugsSameStartEndGpOne$Freq)
	  colnames(dat) <- c("drug1","drug2","Freq")
	  drugsStartSwitchGpOne <- as.data.frame(table(cohortStartSwitchGpOne$drugPath))
	  drugsStartSwitchGpOne <- drugsStartSwitchGpOne[order(-drugsStartSwitchGpOne$Freq),]
	  dat2 <- data.frame(do.call('rbind', strsplit(as.character(drugsStartSwitchGpOne$Var1),'->',fixed=TRUE)))
	  dat2 <- cbind(dat2,drugsStartSwitchGpOne$Freq)
	  colnames(dat2) <- c("drug1","drug2","Freq")
	  datGpOne <- rbind(dat,dat2)
	  write.table(datGpOne,file="datGpOne") 
	  #Group-2
	  drugsSameStartEndGpTwo <- as.data.frame(table(cohortSameStartEndGpTwo$drugPath))
	  drugsSameStartEndGpTwo <- drugsSameStartEndGpTwo[order(-drugsSameStartEndGpTwo$Freq),]
	  dat <- data.frame(do.call('rbind', strsplit(as.character(drugsSameStartEndGpTwo$Var1),'->',fixed=TRUE)))
	  dat <- cbind(dat,drugsSameStartEndGpTwo$Freq)
	  colnames(dat) <- c("drug1","drug2","Freq")
	  drugsStartSwitchGpTwo <- as.data.frame(table(cohortStartSwitchGpTwo$drugPath))
	  drugsStartSwitchGpTwo <- drugsStartSwitchGpTwo[order(-drugsStartSwitchGpTwo$Freq),]
	  dat2 <- data.frame(do.call('rbind', strsplit(as.character(drugsStartSwitchGpTwo$Var1),'->',fixed=TRUE)))
	  dat2 <- cbind(dat2,drugsStartSwitchGpTwo$Freq)
	  colnames(dat2) <- c("drug1","drug2","Freq")
	  datGpTwo <- rbind(dat,dat2)
	  write.table(datGpTwo,file="datGpTwo") 
	  #Group-3
	  drugsSameStartEndGpThree <- as.data.frame(table(cohortSameStartEndGpThree$drugPath))
	  drugsSameStartEndGpThree <- drugsSameStartEndGpThree[order(-drugsSameStartEndGpThree$Freq),]
	  dat <- data.frame(do.call('rbind', strsplit(as.character(drugsSameStartEndGpThree$Var1),'->',fixed=TRUE)))
	  dat <- cbind(dat,drugsSameStartEndGpThree$Freq)
	  colnames(dat) <- c("drug1","drug2","Freq")
	  drugsStartSwitchGpThree <- as.data.frame(table(cohortStartSwitchGpThree$drugPath))
	  drugsStartSwitchGpThree <- drugsStartSwitchGpThree[order(-drugsStartSwitchGpThree$Freq),]
	  dat2 <- data.frame(do.call('rbind', strsplit(as.character(drugsStartSwitchGpThree$Var1),'->',fixed=TRUE)))
	  dat2 <- cbind(dat2,drugsStartSwitchGpThree$Freq)
	  colnames(dat2) <- c("drug1","drug2","Freq")
	  datGpThree <- rbind(dat,dat2)
	  write.table(datGpThree,file="datGpThree") 
}
#'This function generates keyword and ignore lists based on the expansion of
#'concepts.
#'
#'@description Given any given concept_id or string of text this function
#'generates keyword and ignore lists based on the expansion of concepts (looking
#'at their synonyms).
#'
#'@param connection    The connection to the database server.
#'@param aphroditeConceptName  The string of text / concept name to use.
#'@param schema        The database schema being used.
#'@param dbms          The target DBMS for SQL to be rendered in.
#'
#'@details Takes the aphroditeConceptName looks for synonyms and builds a list
#'of related concepts using the vocabulary hierarchies
#'
#'@return A list with two elements: a list of positive keywords found
#'(keywordlist_ALL), and a list of ignore keywords (ignorelist_ALL)
#'
#' @examples \dontrun{
#'
#'wordLists <- buildKeywordList(conn, aphrodite_concept_name, cdmSchema, dbms)
#'
#' }
#'
#'@export
fireItUp <- function(){
### I am leaving this here for now, but will probably remove it ALL the library stuff should be out in the main call
	library(SqlRender)
	library(data.table)
	library(DatabaseConnector)
	library(plyr)
	library(dplyr)
	library(reshape2)
	#library(randomForest)
	library(lattice)
	library(ggplot2)
	library(caret) # Could not install this on Dev-3. It says package ‘caret’ is not available (for R version 3.1.2)
	library(pROC) # Could not install this on Dev-3. It Says package ‘pROC’ is not available (for R version 3.1.2)
	library(iterators)
	library(parallel)
	library(doMC)
	library(gridExtra)
	#library(glmnet)
	library(splines)
	library(survival)
	library(Rmisc)

	#This calls all extra stuff
	source('~/R/aux_functions.R')

	cdmSchema <- "ohdsiv5"
	resultsSchema <- "results_schema"
	sourceName <- "source_name"
	dbms <- "postgresql" #Should be "sql server", "oracle", "postgresql" or "redshift"

	user <- "achilles" #Change when in production to: NULL
	pw <- "achilles" #Change when in production to: NULL
	server <- "localhost/jmbanda"  #Change when in production to: "server_name"
	port <- "5432" #Change when in production to: NULL


	#Create database connection - this should take values from a settings file for easy configuration - see aphrodite code
	jdbcDrivers <<- new.env()
	connectionDetails <- createConnectionDetails(dbms="postgresql", server="localhost/jmbanda", user="achilles", password="achilles" ,schema="ohdsiv5", port="5432")
	conn <- connect(connectionDetails)

	#------------------------------ Building the Cohort --------------------------------------
	message("-------- The Overall Process May take an Hour, But not more than that ---------
			--------------------------------------------------------------------------------")

	message("********** Building Cohort **********")
	source('~/OHDSI-Final/OHDSI-FinalScripts/getCohort.R')
	getCohort(conn)
	# This is one major script that will query the database. The cohort output is saved in the current directory, which will be used for
	# most of the subsequent calculations

	#------------------------------ Getting Basic Stats --------------------------------------
	message("********* Getting Basic Statistics *********")
	source('~/OHDSI-Final/OHDSI-FinalScripts/getBasicStatistics.R')
	getBasicStatistics()

	#----------------------------- Getting Treatment Pathways -------------------------------
	message("********* Building Treatment Pathways ********")
	source('~/OHDSI-Final/OHDSI-FinalScripts/getTreatmentPath.R')
	getTreatmentPath()

	#---------------------------- Getting HbA1c Profile Distribution ------------------------
	message("********* Getting HbA1c Distribution After Treatment **********")
	source('~/OHDSI-Final/OHDSI-FinalScripts/getHbADisb.R')
	getHbADisb()

	#---------------------------- Learning Features Before Treatment ------------------------
	# This function too uses database query
	message("******** Learning Potential Features Pre-Treatment **********")
	source('~/OHDSI-Final/OHDSI-FinalScripts/getCombinedFeatures.R')
	getCombinedFeatures(conn)

	#--------------------------- Perform Survival Analysis ---------------------------------
	message("********* Performing Survival Analysis ************")
	source('~/OHDSI-Final/OHDSI-FinalScripts/getSurvival.R')
	getSurvival()
}