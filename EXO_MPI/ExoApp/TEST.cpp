int isAllFinish=0, isMySendFinish=0;

for (int i=0; i<m; i++>) {
  /* Envoi synchrone */
  MPI_Issend(sbuf[i],size[i],type,dst[i],tag,comm,&reqs[i]);
 }
while (!isAllFinish) {
  /* Avons nous un message pret a etre recu */
  MPI_Iprobe(MPI_ANY_SOURCE,tag,comm,&flag,&astat);
  if (flag) {
    /* Recoit le message */
    MPI_Recv(rbuf,rsize,type,astat%MPI_SOURCE,tag,comm,&rstat);
  }
  if (!isMySendFinish) {
    /* Verifie si tous nos ssend sont termines */
    MPI_Testall(m,reqs,&flag,MPI_STATUSES_IGNORE);
    if (flag) {
      /* Si c’est le cas on demarre la ibarrier */
      MPI_Ibarrier(comm,&reqb);
      isMySendFinish=1;
    }
  } else {
    /* Test si tout le monde a fait la ibarrier */
    MPI_Test(&reqb,isAllFinish,MPI_STATUS_IGNORE);
  }
 }
