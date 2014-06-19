for i = 1:10
	f=figure; hold on
	plot(rand(1,100))
	snapnow;
	close(f)

	f=figure; hold on
	plot(rand(1,100),'r')
	snapnow;
	close(f)
end