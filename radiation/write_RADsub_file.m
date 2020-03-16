function write_RADsub_file(out_path,tid,SW,LW,infoRAD)
    filename = [out_path,datestr(tid(1),'yyyymm')];
    save(filename,'tid','SW','LW','infoRAD')
end