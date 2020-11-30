function [dart] = GTwave_DART_out(dart,modspace,kappa)

dart.win_length

 wave_out=modspace.Tgrn*modspace.ds;
 for aa=1:length(dart.id)
     FNAME= strcat('Dart',num2str(dart.id(aa)),'_kp',num2str(kappa),'.txt');
     if aa==1
         wave=wave_out(1:dart.win_length(aa,2));
         wave_t=[dart.window(aa,1):dart.win_t(aa,3):dart.window(aa,2)]';
         wtt=[wave,wave_t];
         save(FNAME,'wtt','-ascii')
     else
         wave=wave_out(dart.win_length(aa-1,2):dart.win_length(aa,2));
         wave_t=[dart.window(aa,1):dart.win_t(aa,3):dart.window(aa,2)]';
         wtt=[wave,wave_t];
         save(FNAME,'wtt','-ascii')
     end
 end

end

