clc;
clear all;
close all;


% default values %

interval = 0.5;

t = 0:0.01:5;

heartRate = 72;
pwav  = [0.25 0.09 0.16];
qwav  = [0.025 0.066 0.166];
qrswav  = [1.6 0.11];
swav  = [0.25 0.066 0.09];
twav  = [0.35 0.142 0.2];
uwav  = [0.035 0.0476 0.433];

% default values %

default = input('Default values for ecg signal? [y/n] : ', 's'); 

if(default ~= 'y')

  heartRate = input('Enter the heart beat rate : ');

  %p wave specifications
  fprintf('\n\np wave specifications\n');
  d = input('Default values ? [y/n] : ', 's');
  if(d ~= 'y')
     pwav(1) = input('amplitude = ');
     pwav(2) = input('duration = ');
     pwav(3) = input('p-r interval = ');
  end


  %q wave specifications
  fprintf('\n\nq wave specifications\n');
  d = input('Default values ? [y/n] : ', 's'); 
  if(d ~= 'y')
     qwav(1) = input('amplitude = ');
     qwav(2) = input('duration = ');
     qwav(3) = 0.166;
  end


  %qrs wave specifications
  fprintf('\n\nqrs wave specifications\n');
  d = input('Default values ? [y/n] : ', 's');
  if(d ~= 'y')
     qrswav(1) = input('amplitude = ');
     qrswav(2) = input('duration = ');
  end



  %s wave specifications
  fprintf('\n\ns wave specifications\n');
  d = input('Default values ? [y/n] : ', 's'); 
  if(d ~= 'y')
     swav(1) = input('amplitude = ');
     swav(2) = input('duration = ');
     swav(3) = 0.09;
  end


  %t wave specifications
  fprintf('\n\nt wave specifications\n');
  d = input('Default values ? [y/n] : ', 's'); 
  if(d ~= 'y')
     twav(1) = input('amplitude = ');
     twav(2) = input('duration = ');
     twav(3) = input('s-t interval = ');
  end

  
  %u wave specifications
  fprintf('\n\nu wave specifications\n');
  d=input('Default values ? [y/n] : ', 's'); 
  if(d ~= 'y')
     uwav(1) = input('amplitude = ');
     uwav(2) = input('duration = ');
     uwav(3) = 0.433;
  end
end


li = 30 / heartRate;

fprintf('\n\nECG Simulator is running...');
fprintf('\n[ Press Contol + C to stop ]\n');

figure('name', 'ECG Simulation', 'numbertitle', 'off','Units','normalized','Position',[0.02 0.4 0.9 0.5]);

while(true)
  
  t = t+ 0.25;
  li=rand*li;
  %pwav output
  pwav_result = p_wav(t, pwav(1), pwav(2), pwav(3), li);

  %qwav output
  qwav_result = q_wav(t, qwav(1), qwav(2), qwav(3), li);

  %qrswav output
  qrswav_result = qrs_wav(t, qrswav(1), qrswav(2), li);

  %swav output
  swav_result = s_wav(t, swav(1), swav(2), swav(3), li);

  %twav output
  twav_result = t_wav(t, twav(1), twav(2), twav(3), li);

  %uwav output
  uwav_result = u_wav(t, uwav(1), uwav(2), uwav(3), li);

  %ecg output
  ecg = pwav_result + qwav_result + qrswav_result;
  ecg =ecg+ swav_result + twav_result + uwav_result;
  

  plot(t, ecg);
  grid minor;
  xlim([min(t), max(t)]);
  
  pause(interval);
end
