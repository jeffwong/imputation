projection.time = function(timestamps) {
  year = year(timestamps)
  month = month(timestamps)
  yday = yday(timestamps)
  mday = mday(timestamps)
  hour = hour(timestamps)
  minute = minute(timestamps)
  weekday = wday(timestamps)
  weekend = (weekday == 1 | weekday == 7)
  list(year = year, month = month, yday = yday,
       mday = mday, hour = hour, minute = minute,
       weekday = weekday, weekend = weekend)
}
