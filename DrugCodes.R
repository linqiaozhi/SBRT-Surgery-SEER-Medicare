#Drug Codes of Interest (note that chemotherapy pertains to chemotherapies typically given to lung cancer patients)

#Chemotherapies 
chemotherapy_all.dxs<-c( )
carboplatin.dxs<-c()
cisplatin.dxs<-c()


#Anticoagulants 
anticoagulants_all.dxs<-c()
heparin.dxs<-c()
rivaroxban

#Insulins 



#Add 0s to the front (example code)

novo.complete.dxs<-list()
for (i in 1:length(novo.dxs))
{novo.complete.dxs[i]<-ifelse(nchar(novo.dxs[i])<11, paste('00', novo.dxs[i], sep=""), novo.dxs[i])

}

novo.complete.dxs


carbo.complete.dxs