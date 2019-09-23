function h = light_on(t, pulse_duration, n_pulses, pulse_interval, train_interval)

train_duration = n_pulses*(pulse_duration + pulse_interval) + train_interval;
t = rem(t,train_duration);
if t > n_pulses*(pulse_duration + pulse_interval)
    h = 0;
    return
elseif rem(t,pulse_duration+pulse_interval)>pulse_duration
    h = 0;
else
    h = 1;
    return
end
   