void make_script()
{
  string name = "/SingleNeutrino_Pt-2To20-gun/SingleNeutrino_Pt-2To20_Run3Winter21DRMiniAOD_01092021/210901_092535/000";
  int nfile = 7052;
  cout<<"ssd "<<name<<nfile<<endl;
  ofstream f;
  f.open("copy.sh");
  for(int j=1;j<=nfile;j++)
    {
      int t=j/1000;
      string temp="xrdcp root://se01.indiacms.res.in//store/user/rsaxena"+name+to_string(t)+"/RoorTree_reEmulL1_MC_RAW_RelValZTT_20210727_"+to_string(j)+".root .";
      f<<temp<<endl;
    }
  f.close();
}
