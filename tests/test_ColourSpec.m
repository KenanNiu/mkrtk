function test_ColourSpec


rgb = [...
    1 0 0
    0 1 0
    0 0 1
    1 1 0
    0 1 1
    1 0 1
    0 0 0
    1 1 1];

snames = {'r','g','b','y','c','m','k','w'};
lnames = {'red','green','blue','yellow','cyan','magenta','black','white'};
mixed = {'r',[0 1 0],'blue','y',[0 1 1],'magenta','k',[1 1 1]};
assertEqual(rgb,ColourSpec('torgb',lnames))
assertEqual(lnames,ColourSpec('tolongnames',mixed))
assertEqual(lnames,ColourSpec('tolongnames',rgb))
assertEqual(lnames,ColourSpec('tolongnames',snames))


