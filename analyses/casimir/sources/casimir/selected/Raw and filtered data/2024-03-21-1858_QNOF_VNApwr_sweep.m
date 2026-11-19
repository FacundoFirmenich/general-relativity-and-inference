%resetgpib;
clear all;
close all;
cd('C:/NEMSmeasure');

tic;
filename=todayfilename;
filename=strcat(filename(1:4),'-',filename(5:6),'-',filename(7:8),'-',filename(9:12))
path=createdir;


OpenGPIB;

% copies the working file
exec_file_name = mfilename;
copyfile([exec_file_name, '.m'], [path,'\',filename,'_',exec_file_name, '.m']);

common.clock=clock;
vari = ['r', 'g', 'b', 'c', 'm', 'y', 'k'];
SpecAn  = sa;
MWsource = SMA100B;
VNA = Instruments.VNA_ZNB8(9);
fprintf(VNA.hdl,'SWE:TIME 1');
VNA.set_output_on();

VNApowers  = [-30:2.5:-10]; %Powers of the VNA drive
MWpowers   = [-20:2.5:-15]; %Powers of the microwave drive
SMA100BSetPower(MWsource, -20)
fcavity = 5.461831e9; %Cavity frequency
f_m = 10.011e6; %Mechanical frequency

% this is how long fft analyzer takes data in single sweep. The longer the
% faster usually. Normal 1s...10s. It matters in readout routine pause times.
sweeptime = 2;
RBW = 30;
VNA_BW = 100;
VNA_dwell_time = 0.05;
spaspan = 100e3;
VNAspan = 30e3;
numfreqs  = 501;



Navg = 1;

SPA2692ASignalAnalyzerMode(SpecAn, fcavity - f_m ,RBW, spaspan, sweeptime);
pause(0.2); % this freaks out sometimes
avgtime = 12; %Time of measurement. Takes averages for %avgtime% seconds

VNA.set_sweep_CW();
VNA.set_BW(VNA_BW);
VNA.set_power(VNApowers(1));


% Determine the frequency points for the VNA
allfreqs_red  = linspace(fcavity-f_m+VNAspan/2,fcavity-f_m-VNAspan/2,numfreqs);
allfreqs_blue = linspace(fcavity+f_m+VNAspan/2,fcavity+f_m-VNAspan/2,numfreqs);


VNA.segment_delete();
for findex=1:numfreqs
    VNA.segment_add(allfreqs_red(findex), allfreqs_red(findex), 1, VNApowers(1), VNA_dwell_time, VNA_BW)
end

VNA.set_sweep_SEGM();

numpoints=str2double(query(SpecAn,'SWE:POIN?'));

%Get the VNA sweep time
VNAsweeptime = VNA_dwell_time*numfreqs;

%Set the SA sweep time to be the same as the VNA sweep time

fprintf(SpecAn, ['CALC:ATIM:LENG ',num2str(VNAsweeptime*1.3) 's']);




frequency_centers = [fcavity-f_m, fcavity+f_m, fcavity-2*f_m,  fcavity+2*f_m, fcavity-3*f_m,  fcavity+3*f_m];

names = ["red sideband", "blue sideband", "second red sideband", "second blue sideband", "third red sideband", "third blue sideband"];

data_up   = zeros([numpoints,length(frequency_centers),length(VNApowers), length(MWpowers), Navg]); % VNA frequency forward sweep
data_down = zeros([numpoints,length(frequency_centers),length(VNApowers), length(MWpowers), Navg]); % VNA frequency backward sweep

%data_red1 = zeros([numpoints,length(VNApowers)]);
%data_red2 = zeros([numpoints,length(VNApowers)]);
%data_red3 = zeros([numpoints,length(VNApowers)]);

%data_blue1 = zeros([numpoints,length(VNApowers)]);
%data_blue2 = zeros([numpoints,length(VNApowers)]);
%data_blue3 = zeros([numpoints,length(VNApowers)]);

freq = zeros([numpoints,length(frequency_centers)]);

for i = 1:length(frequency_centers)
    freq(:,i)= linspace(frequency_centers(i)-spaspan/2,frequency_centers(i)+spaspan/2,numpoints);
end


%freqstore = [linspace(fcavity-f_m-spaspan/2,fcavity-f_m+spaspan/2,numpoints),
%    linspace(fcavity-2*f_m-spaspan/2,fcavity-2*f_m+spaspan/2,numpoints),
%    linspace(fcavity-3*f_m-spaspan/2,fcavity-3*f_m+spaspan/2,numpoints),
%    linspace(fcavity+f_m-spaspan/2,fcavity+f_m+spaspan/2,numpoints),
%    linspace(fcavity+2*f_m-spaspan/2,fcavity+2*f_m+spaspan/2,numpoints),
%    linspace(fcavity+3*f_m-spaspan/2,fcavity+3*f_m+spaspan/2,numpoints)];

% Set GPIB time-out long enough to wait for entire data (some reasonable
% time, because Ctrl+C doesnt work while waiting for GPIB).
set(SpecAn,'Timeout',60);



for iter3 = 1:Navg
    for iter2 = 1:length(MWpowers)
        SMA100BSetPower(MWsource, MWpowers(iter2));
        for iter1 = 1:length(VNApowers)
            VNApower = VNApowers(iter1);
            VNA.set_power(VNApower)
              
            for icenter = 1:length(frequency_centers)
                 SPA2692AsetCenterFreq(SpecAn, frequency_centers(icenter) ,RBW, spaspan);
                 pause(3);
            
                %% Measure from low to high freq
                %Create frequency segments for the VNA so it can sweep up
                VNA.segment_delete()
                for findex=1:numfreqs
                    VNA.segment_add(allfreqs_red(findex), allfreqs_red(findex), 1, VNApower, VNA_dwell_time, VNA_BW)
                end
               
                VNA.set_sweep_SEGM()
            
                fprintf(VNA.hdl,'INIT:CONT OFF');
                fprintf(VNA.hdl,'AVER:CLE'); % reset the average
            
                %Initialize the VNA and spectrum analyzer immediately afterwards
                fprintf(VNA.hdl,'INIT:IMM');
                fprintf(SpecAn,'INIT:MODE:SING');
                
            
                
                % Wait for trace to finish (note: only waits when executing the next command)
                fprintf(SpecAn,'*WAI');
            
                % Read in batches (5121 points maximum)
                data_temp = nan(numpoints,1);
                for dindex = 1:5121:numpoints
                    damount = min(5121, numpoints - dindex + 1);
                    % TRAC? start, stop
                    fprintf(sa,sprintf('TRAC? %d, %d', dindex-1, damount));
                    data_temp(dindex:(dindex+damount-1),1) = scanstr(sa,',','%f');
                end
            
                data_up(:, icenter, iter1,iter2,iter3) = data_temp;
        
        
                %% Measure from high to low freq
                 %Create frequency segments for the VNA so it can sweep up
                VNA.segment_delete()
        
                for findex=numfreqs:-1:1
                    VNA.segment_add(allfreqs_red(findex), allfreqs_red(findex), 1, VNApower, VNA_dwell_time, VNA_BW)
                end
                VNA.set_sweep_SEGM()
            
                fprintf(VNA.hdl,'INIT:CONT OFF');
                fprintf(VNA.hdl,'AVER:CLE'); % reset the average
            
                %Initialize the VNA and spectrum analyzer immediately afterwards
                fprintf(VNA.hdl,'INIT:IMM');
                fprintf(SpecAn,'INIT:MODE:SING');
                
            
                
                % Wait for trace to finish (note: only waits when executing the next command)
                fprintf(SpecAn,'*WAI');
            
                % Read in batches (5121 points maximum)
                data_temp = nan(numpoints,1);
                for dindex = 1:5121:numpoints
                    damount = min(5121, numpoints - dindex + 1);
                    % TRAC? start, stop
                    fprintf(sa,sprintf('TRAC? %d, %d', dindex-1, damount));
                    data_temp(dindex:(dindex+damount-1),1) = scanstr(sa,',','%f');
                end
            
                data_down(:, icenter, iter1,iter2,iter3) = data_temp;
        
        
            %% plotting
        
            try
                figure(1)
                subplot(length(frequency_centers),1, icenter)
                plot(freq(:,icenter), data_up(:, icenter, iter1,iter2,iter3))
                hold on;
                plot(fliplr(freq(:,icenter)), data_down(:, icenter, iter1,iter2,iter3))
                title(names(icenter));
            catch e
            end
          
        
            end
        
        
        
        end
    
    
    end
end







save([path,'\',filename]);
fprintf(MWsource, 'OUTP OFF');
