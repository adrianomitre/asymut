#!/usr/bin/env ruby

$text = STDIN.readlines

def compare(a, b)
  return 1 if a =~ /End_of_file/
  return -1 if b =~ /End_of_file/

  a_ch, a_time = a.split(',')[0..1].collect {|val| val.to_i }
  b_ch, b_time = b.split(',')[0..1].collect {|val| val.to_i }
  if a_ch != b_ch
    a_ch.to_i <=> b_ch
  elsif a_time != b_time
    a_time <=> b_time
  else
    $text.index(a) <=> $text.index(b)
  end
end

puts $text.sort {|a, b| compare(a, b) }
