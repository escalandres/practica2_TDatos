function crcbit = FCRC(H)

bS = "";
for i=1:length(H)
    a = string(H(i));
    %bS = horzcat(bS,a);
    bS = bS+a;
end

smessage="0b"+bS+"u32";
%message = 0b1101100111011010u32;
message = str2num(smessage);
messageLength = length(H);
divisor = 0b1101u32;
divisorDegree= 3;

divisor = bitshift(divisor,messageLength-divisorDegree-1);
dec2bin(divisor)

divisor = bitshift(divisor,divisorDegree);
remainder = bitshift(message,divisorDegree);
dec2bin(divisor)

dec2bin(remainder)

for k=1:messageLength
	if bitget(remainder,messageLength+divisorDegree)
		remainder = bitxor(remainder,divisor);
	end
	remainder = bitshift(remainder,1);
end

CRC_check_value = bitshift(remainder,-messageLength);
FCS = dec2bin(CRC_check_value);
H(length(H)+1)=FCS(1);
H(length(H)+1)=FCS(2);
H(length(H)+1)=FCS(3);
crcbit = H;

end