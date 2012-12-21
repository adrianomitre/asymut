#!/usr/bin/env ruby
$stdin.each do |line|
  fields = line.split
  if fields.size <= 3 # mono
    puts line
  elsif fields[1].to_i == 1
    puts fields.values_at(0,2,3).join(" ") # discards channel number
  end
end
