file1 = filepath('Night.jpg', $
  subdirectory=['examples', 'data'])

im1 = image(file1, $
  background_color='midnight blue', $
  image_dimensions=[360, 180], $
  image_location=[-180, -90], $
  xrange=[-180, 0], $
  yrange=[-90, 90], $
  dimensions=[512, 512], $ ; dimensions for plotting image area
  margin=0)

file2 = filepath('Day.jpg',  $
  subdirectory=['examples', 'data'])

im2 = image(file2, $
  /overplot, $
  image_dimensions=[360, 180], $
  image_location=[-180, -90], $
  transparency=50)

t = text(-175, 80, $
  '$it Day/Night$', $
  /data, $
  font_size=20)

for i=-100, 100 do im2.transparency=abs(i)
end




