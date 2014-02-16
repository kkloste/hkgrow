function recipient = matlabmail(recipient, message, subject, sender, psswd)
%MATLABMAIL( recipient, message, sender, psswd)
%   sends the character string stored in 'message' with subjectline 'subject'
%   to the address in 'recipient', from the email address 'sender'.
%   This requires that the sending address is a GMAIL email account.
%
%   Note: Authentication failed when my gmail account had
%   2-step verification enabled.
%
% Example:
% matlabmail('kyle.kloster@gmail.com', 'experiment [name] done', 'matlab - [xp name] done');



if nargin<5
    sender = 'puma.notifs@gmail.com';
    psswd = 'JrAAmeBFUpza';
end

setpref('Internet','E_mail',sender);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',sender);
setpref('Internet','SMTP_Password',psswd);

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', ...
                  'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

sendmail(recipient, subject, message);
end

