%  This script constructs file names list (returned in fileNamesMT)
%		for one site
%  It is called after ZplotWhat, which is edited
%  to control what is plotted (site, band, time window, etc.)
%   All variables required are defined in ZplotWhat ... names
%    of these variables have to agree with what is in this script 

%  Local and remote file names
LocSite = SITES{iSite};
RemSite = SITES{3-iSite};

if RR == 1
   subDir = 'MTR';
   suffix = [SUFFIX '.zrr'];
   root = ['RR' SUFFIX];
elseif RR == 0
   subDir = 'MTS';
   suffix = [SUFFIX '.zss'];
   root = ['SS' SUFFIX];
else
   subDir = 'MMT';
   suffix = [SUFFIX '.zmm'];
   root = ['MM' SUFFIX];
end

if BAND == 40
   %  make list of starting hours
   hr = 0;
   for k = 1:nHour
     H = num2str(hr);
     if(hr < 10)
       H = ['0' H];
     end
     HRS{k} = H;
     hr = hr+hourStep;
   end
end

%  compute last day, numbered consecutively from yeardays in
%   first year
lastDay = dayEnd;
for yr = YearEnd:-1:YearOne+1
   if(mod(yr,4) == 0 )
      lastDay = lastDay+366;
   else
      lastDay = lastDay+365;
   end
end

k = 0;
yr = YearOne;
nDaySub = 0;
if(mod(yr,4) == 0)
  nDayYr = 366;
else
  nDayYr = 365;
end

clear fileNamesMT;

for DAY = dayOne:dayStep:lastDay
   day = DAY-nDaySub;
   if day > nDayYr
      day = day-nDayYr;
      yr = yr + 1;
      nDaySub = nDaySub+nDayYr;
      if(mod(yr,4) == 0)
        nDayYr = 366;
      else
        nDayYr = 365;
      end
   end
   if(day<10)
      Day = ['00' num2str(day)];  
   elseif (day<100)
      Day = ['0' num2str(day)];  
   else
      Day = num2str(day);  
   end

   if ADDyear
      Day = [Day '_' num2str(yr) ];
   end

   if BAND == 1
      k = k+1;
      type = suffix(end-3:end);
      if(type == '.zrr')
         siteArray =  [LocSite '_' Day  '_' RemSite];
      elseif(type == '.zmm')
         siteArray = [LocSite '_' Day ];
      elseif(type == '.zss')
         siteArray = [LocSite '_' Day ];
      elseif(type == '.stf')
         siteArray = [LocSite '_' Day ];
      end
      fileNamesMT{k} =  ...
         [ PKDSAOdata num2str(yr) '/' subDir '/' siteArray suffix];
    else
       for h = 1:nHour
          k = k+1;
          type = suffix(end-3:end);
          if(type == '.zrr')
             siteArray =  [LocSite '_' Day  '_' HRS{h} '_' RemSite];
          elseif(type == '.zmm')
             siteArray = [LocSite '_' Day '_' HRS{h}];
          elseif(type == '.zss')
             siteArray = [LocSite '_' Day '_' HRS{h}];
          elseif(type == '.stf')
             siteArray = [LocSite '_' Day '_' HRS{h}];
          end
          fileNamesMT{k} =  ...
             [ PKDSAOdata num2str(yr) '/' subDir '/' siteArray suffix];
      end
   end
end
